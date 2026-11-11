// Extended PDF parameter setup; canonical source is setup_clubb_pdf_params.F90.
// Local storage indices are zero-based. Shared metadata/index_mapping IDs retain
// Native ordinals; adapt only at metadata captures and lookup calls.
#include "scratch_workspace.hpp"
#include "setup_clubb_pdf_params.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "error_code.hpp"
#include "pdf_utilities.hpp"
#include "clip_explicit.hpp"
#include "advance_helper_module.hpp"
#include "matrix_operations.hpp"
#include "Nc_Ncn_eqns.hpp"
#include "precipitation_fraction.hpp"
#include "index_mapping.hpp"
#include "diagnose_correlations_module.hpp"
#include "stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
namespace clubb_kokkos {
namespace setup_clubb_pdf_params_detail {
// Fortran TRIM removes trailing blanks only. This character-language adapter
// preserves caller-provided hydrometeor names and dynamic stats names.
std::string trim(std::string value) {
  value.erase(value.find_last_not_of(' ') + 1);
  return value;
}
} // namespace setup_clubb_pdf_params_detail
namespace {
// The component correlations of w and r_t and the component correlations of
// w and theta_l must both be 0 when using the ADG1 PDF.  In other words, w
// and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'), but the
// individual component covariance and correlation are defined to be 0.
// Since the component covariances (or correlations) of w and chi (old s) and
// of w and eta (old t) are based on the covariances (or correlations) of w
// and r_t and of w and theta_l, the individual component correlation and
// covariance of w and chi, as well of as w and eta, are defined to be 0.
//
// The PDF component correlation of w and x, calculated as part of the PDF
// parameters, comes out to be 0 (within numerical roundoff) when the ADG1
// PDF is used.  However, when l_fix_w_chi_eta_correlations is enabled, the
// component correlations of PDF variables are specified, and the values
// that were calculated as part of the PDF parameters are ignored in the
// correlation array.  When this happens, and when the ADG1 PDF is selected,
// the l_follow_ADG1_PDF_standards flag can be used to force the component
// correlations of w and x to have a value of 0, following the ADG1 standard,
// rather than whatever value might be specified in the correlation array.
//
// Note:  the ADG2 PDF follows the same standards as the ADG1 PDF.
constexpr bool l_follow_ADG1_PDF_standards = true;
void calc_corr_norm_and_cholesky_factor(ScratchWorkspace& workspace,                      // InOut
                                        int nzt, int ngrdcol, int iiPDF_type,             // In
                                        int pdf_dim, const hm_metadata_type& hm_metadata, // In
                                        Real2D rc_1, Real2D rc_2,                         // In
                                        Real2D corr_array_n_cloud, Real2D corr_array_n_below, // In
                                        Real4D corr_array_1_n, Real4D corr_array_2_n,         // Out
                                        Real4D corr_cholesky_mtx_1,
                                        Real4D corr_cholesky_mtx_2 // Out
);
void component_corr_w_x(int nzt, int ngrdcol,                               // In
                        Real2D rc_1,                                        // In
                        Real2D rc_2,                                        // In
                        real_t corr_w_x_NN_cloud, real_t corr_w_x_NN_below, // In
                        int iiPDF_type,                                     // In
                        Real2D corr_w_x_1,                                  // Out
                        Real2D corr_w_x_2                                   // Out
);
void component_corr_chi_eta(int nzt, int ngrdcol,         // In
                            Real2D rc_1,                  // In
                            Real2D rc_2,                  // In
                            real_t corr_chi_eta_NN_cloud, // In
                            real_t corr_chi_eta_NN_below, // In
                            bool l_limit_corr_chi_eta,    // In
                            Real2D corr_chi_eta_1,        // Out
                            Real2D corr_chi_eta_2         // Out
);
void component_corr_w_hm_n_ip(int nzt, int ngrdcol,                 // In
                              Real2D corr_w_hm_1_n_in, Real2D rc_1, // In
                              Real2D corr_w_hm_2_n_in, Real2D rc_2, // In
                              real_t corr_w_hm_n_NL_cloud,          // In
                              real_t corr_w_hm_n_NL_below,          // In
                              bool l_calc_w_corr,                   // In
                              Real2D corr_w_hm_1_n,                 // Out
                              Real2D corr_w_hm_2_n                  // Out
);
void component_corr_x_hm_n_ip(int nzt, int ngrdcol,        // In
                              Real2D rc_1,                 // In
                              Real2D rc_2,                 // In
                              real_t corr_x_hm_n_NL_cloud, // In
                              real_t corr_x_hm_n_NL_below, // In
                              Real2D corr_x_hm_1_n,        // Out
                              Real2D corr_x_hm_2_n         // Out
);
void component_corr_hmx_hmy_n_ip(int nzt, int ngrdcol,           // In
                                 Real2D rc_1,                    // In
                                 Real2D rc_2,                    // In
                                 real_t corr_hmx_hmy_n_LL_cloud, // In
                                 real_t corr_hmx_hmy_n_LL_below, // In
                                 Real2D corr_hmx_hmy_1_n,        // Out
                                 Real2D corr_hmx_hmy_2_n         // Out
);
void component_corr_eta_hm_n_ip(int nzt, int ngrdcol,   // In
                                Real2D corr_chi_eta_1,  // In
                                Real2D corr_chi_hm_n_1, // In
                                Real2D corr_chi_eta_2,  // In
                                Real2D corr_chi_hm_n_2, // In
                                Real2D corr_eta_hm_n_1, // Out
                                Real2D corr_eta_hm_n_2  // Out
);
KOKKOS_INLINE_FUNCTION void calc_corr_w_hm_n(real_t wm, real_t wphydrometp,               // In
                                             real_t mu_w_1, real_t mu_w_2,                // In
                                             real_t mu_hm_1, real_t mu_hm_2,              // In
                                             real_t sigma_w_1, real_t sigma_w_2,          // In
                                             real_t sigma_hm_1, real_t sigma_hm_2,        // In
                                             real_t sigma_hm_1_n, real_t sigma_hm_2_n,    // In
                                             real_t mixt_frac, real_t precip_frac_1,      // In
                                             real_t precip_frac_2, real_t hm_tol,         // In
                                             real_t& corr_w_hm_1_n, real_t& corr_w_hm_2_n // Out
);
void pdf_param_hm_stats(int nzt, int ngrdcol, int pdf_dim, int hydromet_dim, // In
                        const hm_metadata_type& hm_metadata,                 // In
                        Real3D hm_1, Real3D hm_2,                            // In
                        Real3D mu_x_1, Real3D mu_x_2,                        // In
                        Real3D sigma_x_1, Real3D sigma_x_2,                  // In
                        Real4D corr_array_1, Real4D corr_array_2,            // In
                        stats_type& stats                                    // InOut
);
void pdf_param_ln_hm_stats(ScratchWorkspace& workspace, int nzt, int ngrdcol, int pdf_dim,
                           const hm_metadata_type& hm_metadata,    // In
                           Real3D mu_x_1_n, Real3D mu_x_2_n,       // In
                           Real3D sigma_x_1_n, Real3D sigma_x_2_n, // In
                           Real4D corr_array_1_n,                  // In
                           Real4D corr_array_2_n,                  // In
                           stats_type& stats                       // InOut
);
void pack_hydromet_pdf_params(int nzt, int ngrdcol,                                  // In
                              int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                              Real3D hm_1, Real3D hm_2, int pdf_dim, Real3D mu_x_1,  // In
                              Real3D mu_x_2, Real3D sigma_x_1, Real3D sigma_x_2,     // In
                              Real4D corr_array_1, Real4D corr_array_2,              // In
                              hydromet_pdf_parameter_collection& hydromet_pdf_params // Out
);
KOKKOS_INLINE_FUNCTION real_t compute_rtp2_from_chi(real_t sigma_chi_1, real_t sigma_chi_2, // In
                                                    real_t sigma_eta_1, real_t sigma_eta_2, // In
                                                    real_t rt_1, real_t rt_2,               // In
                                                    real_t crt_1, real_t crt_2,             // In
                                                    real_t mixt_frac,                       // In
                                                    real_t corr_chi_eta_1,
                                                    real_t corr_chi_eta_2 // In
);
} // namespace
// Thin private section adapters use the shared array owner and restore outputs.
namespace {
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7>
inline void
calc_corr_norm_and_cholesky_factor(ScratchWorkspace& workspace,                           // InOut
                                   int nzt, int ngrdcol, int iiPDF_type,                  // In
                                   int pdf_dim, const hm_metadata_type& hm_metadata,      // In
                                   Field0 rc_1, Field1 rc_2,                              // In
                                   Field2 corr_array_n_cloud, Field3 corr_array_n_below,  // In
                                   Field4 corr_array_1_n, Field5 corr_array_2_n,          // Out
                                   Field6 corr_cholesky_mtx_1, Field7 corr_cholesky_mtx_2 // Out
) {
  auto workspace_scope                = workspace.scope();
  const auto rc_1_exec                = to_execution_space(workspace, rc_1);
  const auto rc_2_exec                = to_execution_space(workspace, rc_2);
  const auto corr_array_n_cloud_exec  = to_execution_space(workspace, corr_array_n_cloud);
  const auto corr_array_n_below_exec  = to_execution_space(workspace, corr_array_n_below);
  const auto corr_array_1_n_exec      = to_execution_space(workspace, corr_array_1_n, false);
  const auto corr_array_2_n_exec      = to_execution_space(workspace, corr_array_2_n, false);
  const auto corr_cholesky_mtx_1_exec = to_execution_space(workspace, corr_cholesky_mtx_1, false);
  const auto corr_cholesky_mtx_2_exec = to_execution_space(workspace, corr_cholesky_mtx_2, false);
  calc_corr_norm_and_cholesky_factor(
      workspace, nzt, ngrdcol, iiPDF_type, pdf_dim, hm_metadata, rc_1_exec, rc_2_exec,
      corr_array_n_cloud_exec, corr_array_n_below_exec, corr_array_1_n_exec, corr_array_2_n_exec,
      corr_cholesky_mtx_1_exec, corr_cholesky_mtx_2_exec);
  from_execution_space(workspace, corr_array_1_n, corr_array_1_n_exec);
  from_execution_space(workspace, corr_array_2_n, corr_array_2_n_exec);
  from_execution_space(workspace, corr_cholesky_mtx_1, corr_cholesky_mtx_1_exec);
  from_execution_space(workspace, corr_cholesky_mtx_2, corr_cholesky_mtx_2_exec);
}
template <class Field0, class Field1, class Field2, class Field3>
inline void component_corr_w_x(int nzt, int ngrdcol,                               // In
                               Field0 rc_1,                                        // In
                               Field1 rc_2,                                        // In
                               real_t corr_w_x_NN_cloud, real_t corr_w_x_NN_below, // In
                               int iiPDF_type,                                     // In
                               Field2 corr_w_x_1,                                  // Out
                               Field3 corr_w_x_2                                   // Out
) {
  const auto rc_1_exec       = to_execution_space(rc_1);
  const auto rc_2_exec       = to_execution_space(rc_2);
  const auto corr_w_x_1_exec = to_execution_space(corr_w_x_1, false);
  const auto corr_w_x_2_exec = to_execution_space(corr_w_x_2, false);
  component_corr_w_x(nzt, ngrdcol, rc_1_exec, rc_2_exec, corr_w_x_NN_cloud, corr_w_x_NN_below,
                     iiPDF_type, corr_w_x_1_exec, corr_w_x_2_exec);
  from_execution_space(corr_w_x_1, corr_w_x_1_exec);
  from_execution_space(corr_w_x_2, corr_w_x_2_exec);
}
template <class Field0, class Field1, class Field2, class Field3>
inline void component_corr_chi_eta(int nzt, int ngrdcol,         // In
                                   Field0 rc_1,                  // In
                                   Field1 rc_2,                  // In
                                   real_t corr_chi_eta_NN_cloud, // In
                                   real_t corr_chi_eta_NN_below, // In
                                   bool l_limit_corr_chi_eta,    // In
                                   Field2 corr_chi_eta_1,        // Out
                                   Field3 corr_chi_eta_2         // Out
) {
  const auto rc_1_exec           = to_execution_space(rc_1);
  const auto rc_2_exec           = to_execution_space(rc_2);
  const auto corr_chi_eta_1_exec = to_execution_space(corr_chi_eta_1, false);
  const auto corr_chi_eta_2_exec = to_execution_space(corr_chi_eta_2, false);
  component_corr_chi_eta(nzt, ngrdcol, rc_1_exec, rc_2_exec, corr_chi_eta_NN_cloud,
                         corr_chi_eta_NN_below, l_limit_corr_chi_eta, corr_chi_eta_1_exec,
                         corr_chi_eta_2_exec);
  from_execution_space(corr_chi_eta_1, corr_chi_eta_1_exec);
  from_execution_space(corr_chi_eta_2, corr_chi_eta_2_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5>
inline void component_corr_w_hm_n_ip(int nzt, int ngrdcol,                 // In
                                     Field0 corr_w_hm_1_n_in, Field1 rc_1, // In
                                     Field2 corr_w_hm_2_n_in, Field3 rc_2, // In
                                     real_t corr_w_hm_n_NL_cloud,          // In
                                     real_t corr_w_hm_n_NL_below,          // In
                                     bool l_calc_w_corr,                   // In
                                     Field4 corr_w_hm_1_n,                 // Out
                                     Field5 corr_w_hm_2_n                  // Out
) {
  const auto corr_w_hm_1_n_in_exec = to_execution_space(corr_w_hm_1_n_in);
  const auto rc_1_exec             = to_execution_space(rc_1);
  const auto corr_w_hm_2_n_in_exec = to_execution_space(corr_w_hm_2_n_in);
  const auto rc_2_exec             = to_execution_space(rc_2);
  const auto corr_w_hm_1_n_exec    = to_execution_space(corr_w_hm_1_n, false);
  const auto corr_w_hm_2_n_exec    = to_execution_space(corr_w_hm_2_n, false);
  component_corr_w_hm_n_ip(nzt, ngrdcol, corr_w_hm_1_n_in_exec, rc_1_exec, corr_w_hm_2_n_in_exec,
                           rc_2_exec, corr_w_hm_n_NL_cloud, corr_w_hm_n_NL_below, l_calc_w_corr,
                           corr_w_hm_1_n_exec, corr_w_hm_2_n_exec);
  from_execution_space(corr_w_hm_1_n, corr_w_hm_1_n_exec);
  from_execution_space(corr_w_hm_2_n, corr_w_hm_2_n_exec);
}
template <class Field0, class Field1, class Field2, class Field3>
inline void component_corr_x_hm_n_ip(int nzt, int ngrdcol,        // In
                                     Field0 rc_1,                 // In
                                     Field1 rc_2,                 // In
                                     real_t corr_x_hm_n_NL_cloud, // In
                                     real_t corr_x_hm_n_NL_below, // In
                                     Field2 corr_x_hm_1_n,        // Out
                                     Field3 corr_x_hm_2_n         // Out
) {
  const auto rc_1_exec          = to_execution_space(rc_1);
  const auto rc_2_exec          = to_execution_space(rc_2);
  const auto corr_x_hm_1_n_exec = to_execution_space(corr_x_hm_1_n, false);
  const auto corr_x_hm_2_n_exec = to_execution_space(corr_x_hm_2_n, false);
  component_corr_x_hm_n_ip(nzt, ngrdcol, rc_1_exec, rc_2_exec, corr_x_hm_n_NL_cloud,
                           corr_x_hm_n_NL_below, corr_x_hm_1_n_exec, corr_x_hm_2_n_exec);
  from_execution_space(corr_x_hm_1_n, corr_x_hm_1_n_exec);
  from_execution_space(corr_x_hm_2_n, corr_x_hm_2_n_exec);
}
template <class Field0, class Field1, class Field2, class Field3>
inline void component_corr_hmx_hmy_n_ip(int nzt, int ngrdcol,           // In
                                        Field0 rc_1,                    // In
                                        Field1 rc_2,                    // In
                                        real_t corr_hmx_hmy_n_LL_cloud, // In
                                        real_t corr_hmx_hmy_n_LL_below, // In
                                        Field2 corr_hmx_hmy_1_n,        // Out
                                        Field3 corr_hmx_hmy_2_n         // Out
) {
  const auto rc_1_exec             = to_execution_space(rc_1);
  const auto rc_2_exec             = to_execution_space(rc_2);
  const auto corr_hmx_hmy_1_n_exec = to_execution_space(corr_hmx_hmy_1_n, false);
  const auto corr_hmx_hmy_2_n_exec = to_execution_space(corr_hmx_hmy_2_n, false);
  component_corr_hmx_hmy_n_ip(nzt, ngrdcol, rc_1_exec, rc_2_exec, corr_hmx_hmy_n_LL_cloud,
                              corr_hmx_hmy_n_LL_below, corr_hmx_hmy_1_n_exec,
                              corr_hmx_hmy_2_n_exec);
  from_execution_space(corr_hmx_hmy_1_n, corr_hmx_hmy_1_n_exec);
  from_execution_space(corr_hmx_hmy_2_n, corr_hmx_hmy_2_n_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5>
inline void component_corr_eta_hm_n_ip(int nzt, int ngrdcol,   // In
                                       Field0 corr_chi_eta_1,  // In
                                       Field1 corr_chi_hm_n_1, // In
                                       Field2 corr_chi_eta_2,  // In
                                       Field3 corr_chi_hm_n_2, // In
                                       Field4 corr_eta_hm_n_1, // Out
                                       Field5 corr_eta_hm_n_2  // Out
) {
  const auto corr_chi_eta_1_exec  = to_execution_space(corr_chi_eta_1);
  const auto corr_chi_hm_n_1_exec = to_execution_space(corr_chi_hm_n_1);
  const auto corr_chi_eta_2_exec  = to_execution_space(corr_chi_eta_2);
  const auto corr_chi_hm_n_2_exec = to_execution_space(corr_chi_hm_n_2);
  const auto corr_eta_hm_n_1_exec = to_execution_space(corr_eta_hm_n_1, false);
  const auto corr_eta_hm_n_2_exec = to_execution_space(corr_eta_hm_n_2, false);
  component_corr_eta_hm_n_ip(nzt, ngrdcol, corr_chi_eta_1_exec, corr_chi_hm_n_1_exec,
                             corr_chi_eta_2_exec, corr_chi_hm_n_2_exec, corr_eta_hm_n_1_exec,
                             corr_eta_hm_n_2_exec);
  from_execution_space(corr_eta_hm_n_1, corr_eta_hm_n_1_exec);
  from_execution_space(corr_eta_hm_n_2, corr_eta_hm_n_2_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7>
inline void pdf_param_hm_stats(int nzt, int ngrdcol, int pdf_dim, int hydromet_dim, // In
                               const hm_metadata_type& hm_metadata,                 // In
                               Field0 hm_1, Field1 hm_2,                            // In
                               Field2 mu_x_1, Field3 mu_x_2,                        // In
                               Field4 sigma_x_1, Field5 sigma_x_2,                  // In
                               Field6 corr_array_1, Field7 corr_array_2,            // In
                               stats_type& stats                                    // InOut
) {
  const auto hm_1_exec         = to_execution_space(hm_1);
  const auto hm_2_exec         = to_execution_space(hm_2);
  const auto mu_x_1_exec       = to_execution_space(mu_x_1);
  const auto mu_x_2_exec       = to_execution_space(mu_x_2);
  const auto sigma_x_1_exec    = to_execution_space(sigma_x_1);
  const auto sigma_x_2_exec    = to_execution_space(sigma_x_2);
  const auto corr_array_1_exec = to_execution_space(corr_array_1);
  const auto corr_array_2_exec = to_execution_space(corr_array_2);
  pdf_param_hm_stats(nzt, ngrdcol, pdf_dim, hydromet_dim, hm_metadata, hm_1_exec, hm_2_exec,
                     mu_x_1_exec, mu_x_2_exec, sigma_x_1_exec, sigma_x_2_exec, corr_array_1_exec,
                     corr_array_2_exec, stats);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5>
inline void pdf_param_ln_hm_stats(ScratchWorkspace& workspace, int nzt, int ngrdcol, int pdf_dim,
                                  const hm_metadata_type& hm_metadata,    // In
                                  Field0 mu_x_1_n, Field1 mu_x_2_n,       // In
                                  Field2 sigma_x_1_n, Field3 sigma_x_2_n, // In
                                  Field4 corr_array_1_n,                  // In
                                  Field5 corr_array_2_n,                  // In
                                  stats_type& stats                       // InOut
) {
  auto workspace_scope           = workspace.scope();
  const auto mu_x_1_n_exec       = to_execution_space(workspace, mu_x_1_n);
  const auto mu_x_2_n_exec       = to_execution_space(workspace, mu_x_2_n);
  const auto sigma_x_1_n_exec    = to_execution_space(workspace, sigma_x_1_n);
  const auto sigma_x_2_n_exec    = to_execution_space(workspace, sigma_x_2_n);
  const auto corr_array_1_n_exec = to_execution_space(workspace, corr_array_1_n);
  const auto corr_array_2_n_exec = to_execution_space(workspace, corr_array_2_n);
  pdf_param_ln_hm_stats(workspace, nzt, ngrdcol, pdf_dim, hm_metadata, mu_x_1_n_exec, mu_x_2_n_exec,
                        sigma_x_1_n_exec, sigma_x_2_n_exec, corr_array_1_n_exec,
                        corr_array_2_n_exec, stats);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7>
inline void pack_hydromet_pdf_params(int nzt, int ngrdcol,                                  // In
                                     int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                                     Field0 hm_1, Field1 hm_2, int pdf_dim, Field2 mu_x_1,  // In
                                     Field3 mu_x_2, Field4 sigma_x_1, Field5 sigma_x_2,     // In
                                     Field6 corr_array_1, Field7 corr_array_2,              // In
                                     hydromet_pdf_parameter_collection& hydromet_pdf_params // Out
) {
  const auto hm_1_exec         = to_execution_space(hm_1);
  const auto hm_2_exec         = to_execution_space(hm_2);
  const auto mu_x_1_exec       = to_execution_space(mu_x_1);
  const auto mu_x_2_exec       = to_execution_space(mu_x_2);
  const auto sigma_x_1_exec    = to_execution_space(sigma_x_1);
  const auto sigma_x_2_exec    = to_execution_space(sigma_x_2);
  const auto corr_array_1_exec = to_execution_space(corr_array_1);
  const auto corr_array_2_exec = to_execution_space(corr_array_2);
  pack_hydromet_pdf_params(nzt, ngrdcol, hydromet_dim, hm_metadata, hm_1_exec, hm_2_exec, pdf_dim,
                           mu_x_1_exec, mu_x_2_exec, sigma_x_1_exec, sigma_x_2_exec,
                           corr_array_1_exec, corr_array_2_exec, hydromet_pdf_params);
}
} // namespace
// =============================================================================
void setup_pdf_parameters_api(ScratchWorkspace& workspace,                                // InOut
                              const grid& gr, int nzm, int nzt, int ngrdcol, int pdf_dim, // In
                              int hydromet_dim,                                           // In
                              Real2D Nc_in_cloud, Real2D cloud_frac, Real2D Kh_zm,        // In
                              Real2D ice_supersat_frac, Real3D hydromet, Real3D wphydrometp, // In
                              Real2D corr_array_n_cloud, Real2D corr_array_n_below,          // In
                              const hm_metadata_type& hm_metadata,                           // In
                              const pdf_parameter& pdf_params,                               // In
                              Real2D clubb_params,                                           // In
                              int iiPDF_type,                                                // In
                              bool l_use_precip_frac,                                        // In
                              bool l_diagnose_correlations,                                  // In
                              bool l_calc_w_corr,                                            // In
                              bool l_const_Nc_in_cloud,                                      // In
                              bool l_fix_w_chi_eta_correlations,                             // In
                              err_info_type& err_info,                                // InOut
                              Real3D hydrometp2,                                      // Out
                              Real3D mu_x_1_n, Real3D mu_x_2_n,                       // Out
                              Real3D sigma_x_1_n, Real3D sigma_x_2_n,                 // Out
                              Real4D corr_array_1_n, Real4D corr_array_2_n,           // Out
                              Real4D corr_cholesky_mtx_1, Real4D corr_cholesky_mtx_2, // Out
                              precipitation_fractions& precip_fracs,                  // InOut
                              hydromet_pdf_parameter_collection* hydromet_pdf_params, // Out
                              stats_type& stats                                       // InOut
) {
  // Description:
  // Set precipitation fractions, component moments, normal-space correlations
  // and matrix factors for the extended PDF. Pack optional hydrometeor fields
  // and record the Native statistics and ordered diagnostic checks.
  // References:
  //-----------------------------------------------------------------------

  // gr: Grid structure with zero-based physical endpoint slots.
  // nzm: Number of model momentum vertical grid levels
  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // pdf_dim: Number of variables in the correlation array
  // hydromet_dim: Number of hydrometeor species
  // Nc_in_cloud: Mean (in-cloud) cloud droplet conc.       [num/kg]
  // cloud_frac: Cloud fraction                            [-]
  // Kh_zm: Eddy diffusivity coef. on momentum levels [m^2/s]
  // ice_supersat_frac: Ice supersaturation fraction              [-]
  // hydromet: Mean of hydrometeor, hm (overall) (t-levs.) [units]
  // wphydrometp: Covariance < w'h_m' > (momentum levels)     [(m/s)units]
  // corr_array_n_cloud: Prescribed normal space corr. array in cloud  [-]
  // corr_array_n_below: Prescribed normal space corr. array below cl. [-]
  // hm_metadata: Hydrometeor/PDF indices, names, flags, tolerances and variance ratios.
  // pdf_params: PDF parameters                               [units vary]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF.
  // l_use_precip_frac: Flag to use precipitation fraction in KK microphysics. The precipitation
  // fraction is automatically set to 1 when this flag is turned off.
  // l_diagnose_correlations: Diagnose correlations instead of using fixed ones
  // l_calc_w_corr: Calculate the correlations between w and the hydrometeors
  // l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K)
  // l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta)
  // err_info: err_info struct containing err_code and err_header
  // hydrometp2: Variance of a hydrometeor (overall) (m-levs.)   [units^2]
  // mu_x_1_n: Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // mu_x_2_n: Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // corr_cholesky_mtx_1: Transposed corr. cholesky matrix, 1st comp. [-]
  // corr_cholesky_mtx_2: Transposed corr. cholesky matrix, 2nd comp. [-]
  // precip_fracs: Precipitation fractions      [-]
  // hydromet_pdf_params: Hydrometeor PDF parameters        [units vary]
  // stats: Native statistics collector.

  auto workspace_scope = workspace.scope();
  Real2D Kh_zm_c_K_hm =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::Kh_zm_c_K_hm", ngrdcol,
                                   nzm); // Eddy diffusivity coef. on momentum levels [m^2/s]
  Real2D sigma_w_1 =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::sigma_w_1", ngrdcol,
                                   nzt); // Standard deviation of w (1st PDF component)      [m/s]
  Real2D sigma_w_2 =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::sigma_w_2", ngrdcol,
                                   nzt); // Standard deviation of w (2nd PDF component)      [m/s]
  Real2D Ncnm = workspace_scope.view<Real2D>(
      "setup_pdf_parameters_api::Ncnm", ngrdcol,
      nzt); // Mean cloud nuclei concentration, < N_cn >        [num/kg]
  std::string hm_type{};
  std::string var_name{};
  Real2D rcm_pdf = workspace_scope.view<Real2D>("setup_pdf_parameters_api::rcm_pdf", ngrdcol,
                                                nzt); // Liquid water mixing ratio          [kg/kg]
  Real2D wpNcnp_zm = workspace_scope.view<Real2D>(
      "setup_pdf_parameters_api::wpNcnp_zm", ngrdcol,
      nzm); // Covariance of N_cn and w (momentum levs.)   [(m/s)(num/kg)]
  Real2D wpNcnp_zt = workspace_scope.view<Real2D>(
      "setup_pdf_parameters_api::wpNcnp_zt", ngrdcol,
      nzt); // Covariance of N_cn and w on thermo. levels  [(m/s)(num/kg)]
  Real3D hm_1 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::hm_1", ngrdcol, nzt,
      hydromet_dim); // Mean of a precip. hydrometeor (1st PDF component)   [units vary]
  Real3D hm_2 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::hm_2", ngrdcol, nzt,
      hydromet_dim); // Mean of a precip. hydrometeor (2nd PDF component)   [units vary]
  Real3D hydrometp2_zt = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::hydrometp2_zt", ngrdcol, nzt,
      hydromet_dim); // Variance of a hydrometeor (overall); t-lev   [units^2]
  Real3D wphydrometp_zt = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::wphydrometp_zt", ngrdcol, nzt,
      hydromet_dim); // Covariance of w and hm interp. to t-levs. [(m/s)units]
  // Native defines true correlations only when the optional hydro output is present.
  // Their other entries remain outside the defined source contract.
  Real4D corr_array_1 = workspace_scope.view<Real4D>(
      "setup_pdf_parameters_api::corr_array_1", ngrdcol, nzt, pdf_dim,
      pdf_dim); // Correlation array of PDF vars. (comp. 1)             [-]
  Real4D corr_array_2 = workspace_scope.view<Real4D>(
      "setup_pdf_parameters_api::corr_array_2", ngrdcol, nzt, pdf_dim,
      pdf_dim); // Correlation array of PDF vars. (comp. 2)             [-]
  Real3D mu_x_1 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::mu_x_1", ngrdcol, nzt,
      pdf_dim); // Mean array of PDF vars. (1st PDF component)    [units vary]
  Real3D mu_x_2 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::mu_x_2", ngrdcol, nzt,
      pdf_dim); // Mean array of PDF vars. (2nd PDF component)    [units vary]
  Real3D sigma_x_1 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::sigma_x_1", ngrdcol, nzt,
      pdf_dim); // Standard deviation array of PDF vars (comp. 1) [units vary]
  Real3D sigma_x_2 = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::sigma_x_2", ngrdcol, nzt,
      pdf_dim); // Standard deviation array of PDF vars (comp. 2) [units vary]
  Real1D corr_array_scaling =
      workspace_scope.view<Real1D>("setup_pdf_parameters_api::corr_array_scaling", pdf_dim);
  Real3D sigma2_on_mu2_ip_1 =
      workspace_scope.view<Real3D>("setup_pdf_parameters_api::sigma2_on_mu2_ip_1", ngrdcol, nzt,
                                   pdf_dim); // Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
  Real3D sigma2_on_mu2_ip_2 =
      workspace_scope.view<Real3D>("setup_pdf_parameters_api::sigma2_on_mu2_ip_2", ngrdcol, nzt,
                                   pdf_dim); // Ratio array sigma_hm_2^2/mu_hm_2^2             [-]
  Real2D const_Ncnp2_on_Ncnm2 =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::const_Ncnp2_on_Ncnm2", ngrdcol,
                                   nzt); // Prescribed ratio: <Ncn'^2> to <Ncn>^2  [-]
  Real2D stdev_const_Ncnp2_on_Ncnm2 =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::stdev_const_Ncnp2_on_Ncnm2", ngrdcol,
                                   nzt); // Standard deviation of const_Ncnp2_on_Ncnm2
  Real2D const_corr_chi_Ncn =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::const_corr_chi_Ncn", ngrdcol,
                                   nzt); // Prescribed corr. of chi (old s) & Ncn  [-]
  Real2D const_corr_chi_Ncn_n_cloud =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::const_corr_chi_Ncn_n_cloud", ngrdcol,
                                   nzt); // Prescribed corr. of chi & Ncn in cloud [-]
  Real1D precip_frac_tol =
      workspace_scope.view<Real1D>("setup_pdf_parameters_api::precip_frac_tol",
                                   ngrdcol); // Min. precip. frac. when hydromet. present    [-]
  Real3D wphydrometp_chnge = workspace_scope.view<Real3D>(
      "setup_pdf_parameters_api::wphydrometp_chnge", ngrdcol, nzt,
      hydromet_dim); // Change in wphydrometp_zt: covar. clip. [(m/s)units]
  Real2D wm_zt =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::wm_zt", ngrdcol,
                                   nzt); // Mean vertical velocity, <w>, on thermo. levels  [m/s]
  Real2D wp2_zt = workspace_scope.view<Real2D>(
      "setup_pdf_parameters_api::wp2_zt", ngrdcol,
      nzt); // Variance of w, <w'^2> (interp. to t-levs.)      [m^2/s^2]
  Real2D rtp2_zt_from_chi =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::rtp2_zt_from_chi", ngrdcol, nzt);
  Real2D rtp2_zm_from_chi =
      workspace_scope.view<Real2D>("setup_pdf_parameters_api::rtp2_zm_from_chi", ngrdcol, nzm);
  Real1D omicron =
      workspace_scope.view<Real1D>("setup_pdf_parameters_api::omicron",
                                   ngrdcol); // Relative width parameter, omicron = R / Rmax    [-]
  Real1D zeta_vrnce_rat =
      workspace_scope.view<Real1D>("setup_pdf_parameters_api::zeta_vrnce_rat",
                                   ngrdcol); // Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]
  bool l_corr_array_scaling{};
  std::string hydromet_name{}; // Name of a hydrometeor
  int k{};
  int i{};
  int j{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto pdf_params_stdev_chi_1         = pdf_params.stdev_chi_1;
  const auto pdf_params_stdev_chi_2         = pdf_params.stdev_chi_2;
  const auto pdf_params_stdev_eta_1         = pdf_params.stdev_eta_1;
  const auto pdf_params_stdev_eta_2         = pdf_params.stdev_eta_2;
  const auto pdf_params_rt_1                = pdf_params.rt_1;
  const auto pdf_params_rt_2                = pdf_params.rt_2;
  const auto pdf_params_crt_1               = pdf_params.crt_1;
  const auto pdf_params_crt_2               = pdf_params.crt_2;
  const auto pdf_params_mixt_frac           = pdf_params.mixt_frac;
  const auto hm_metadata_iiPDF_chi          = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta          = hm_metadata.iiPDF_eta;
  const auto pdf_params_varnce_w_1          = pdf_params.varnce_w_1;
  const auto pdf_params_varnce_w_2          = pdf_params.varnce_w_2;
  const auto pdf_params_rc_1                = pdf_params.rc_1;
  const auto pdf_params_rc_2                = pdf_params.rc_2;
  const auto pdf_params_cloud_frac_1        = pdf_params.cloud_frac_1;
  const auto hm_metadata_l_mix_rat_hm       = hm_metadata.l_mix_rat_hm;
  const auto hm_metadata_l_frozen_hm        = hm_metadata.l_frozen_hm;
  const auto hm_metadata_hydromet_tol       = hm_metadata.hydromet_tol;
  const auto pdf_params_cloud_frac_2        = pdf_params.cloud_frac_2;
  const auto pdf_params_ice_supersat_frac_1 = pdf_params.ice_supersat_frac_1;
  const auto pdf_params_ice_supersat_frac_2 = pdf_params.ice_supersat_frac_2;
  const auto precip_fracs_precip_frac       = precip_fracs.precip_frac;
  const auto precip_fracs_precip_frac_1     = precip_fracs.precip_frac_1;
  const auto precip_fracs_precip_frac_2     = precip_fracs.precip_frac_2;
  const auto hm_metadata_Ncnp2_on_Ncnm2     = hm_metadata.Ncnp2_on_Ncnm2;
  const auto hm_metadata_iiPDF_Ncn          = hm_metadata.iiPDF_Ncn;
  const auto pdf_params_chi_1               = pdf_params.chi_1;
  const auto pdf_params_chi_2               = pdf_params.chi_2;
  const auto hm_metadata_hmp2_ip_on_hmm2_ip = hm_metadata.hmp2_ip_on_hmm2_ip;
  const auto gr_k_ub_zm                     = gr.k_ub_zm;
  const auto pdf_params_w_1                 = pdf_params.w_1;
  const auto pdf_params_w_2                 = pdf_params.w_2;
  const auto pdf_params_thl_1               = pdf_params.thl_1;
  const auto pdf_params_thl_2               = pdf_params.thl_2;
  const auto hm_metadata_iiPDF_w            = hm_metadata.iiPDF_w;
  const auto hm_metadata_hydromet_list      = hm_metadata.hydromet_list;
  const auto hm_metadata_hydromet_list_host = to_host_space(workspace, hm_metadata_hydromet_list);

  // Check on device, then preserve the first Native species/level/column failure.
  if (clubb_at_least_debug_level_api(0)) {
    Kokkos::View<int***, Kokkos::LayoutLeft> negative =
        workspace_scope.view<Kokkos::View<int***, Kokkos::LayoutLeft>>(
            "setup_pdf_parameters_api::negative", ngrdcol, nzt, hydromet_dim);
    int negative_count = 0;
    Kokkos::parallel_reduce(
        "setup_pdf_parameters_api::negative_check", ((hydromet_dim * ngrdcol * nzt)),
        KOKKOS_LAMBDA(const int cell, int& count) {
          const int i       = ((cell % ngrdcol));
          const int k       = (((cell / ngrdcol) % nzt));
          const int j       = ((cell / (ngrdcol * nzt)));
          negative(i, k, j) = hydromet(i, k, j) < zero_threshold;
          count += negative(i, k, j);
        },
        negative_count);
    if (negative_count != 0) {
      const auto negative_host = to_host_space(workspace, negative);
      const auto hydromet_host = to_host_space(workspace, hydromet);
      const auto names_host    = to_host_space(workspace, hm_metadata.hydromet_list);
      const auto headers_host  = to_host_space(workspace, err_info.err_header);
      for (int j = 0; j < hydromet_dim; ++j)
        for (int k = 0; k < nzt; ++k)
          for (int i = 0; i < ngrdcol; ++i) {
            if (negative_host(i, k, j)) {
              std::string hydromet_name(names_host.data() + ((names_host.extent(0) * j)),
                                        names_host.extent(0));
              hydromet_name.erase(hydromet_name.find_last_not_of(' ') + 1);
              std::ostringstream diagnostic;
              diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              diagnostic.write(headers_host.data() + ((headers_host.extent(0) * i)),
                               headers_host.extent(0));
              diagnostic << '\n'
                         << " at beginning of setup_pdf_parameters_api: " << hydromet_name << " = "
                         << hydromet_host(i, k, j) << " < " << zero_threshold
                         << " at k = " << (((k + 1))) << '\n';
              std::cerr << diagnostic.str();
              const auto codes = err_info.err_code;
              Kokkos::parallel_for(
                  "setup_pdf_parameters_api::first_fatal", 1,
                  KOKKOS_LAMBDA(const int) { codes(i) = clubb_fatal_error; });
              return;
            }
          }
    }
  }
  // Setup some of the PDF parameters
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "setup_pdf_parameters_api::cells_1", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i     = ((cell % extent_1));
          const int k     = ((cell / extent_1));
          sigma_w_1(i, k) = Kokkos::sqrt(pdf_params_varnce_w_1(i, k));
          sigma_w_2(i, k) = Kokkos::sqrt(pdf_params_varnce_w_2(i, k));
        });
  }
  // Compute rcm_pdf for use within SILHS
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "setup_pdf_parameters_api::cells_2", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i   = ((cell % extent_1));
          const int k   = ((cell / extent_1));
          rcm_pdf(i, k) = compute_mean_binormal(pdf_params_rc_1(i, k), pdf_params_rc_2(i, k),
                                                pdf_params_mixt_frac(i, k));
        });
  }
  // Note on hydrometeor PDF shape:
  // To use a single lognormal over the entire grid level, turn off the
  // l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
  // tunable_parameters.in.
  // To use a single delta-lognormal (single lognormal in-precip.), enable the
  // l_use_precip_frac flag and set omicron to 1 and zeta_vrnce_rat to 0 in
  // tunable_parameters.in.
  // Otherwise, with l_use_precip_frac enabled and omicron and zeta_vrnce_rat
  // values that are not 1 and 0, respectively, the PDF shape is a double
  // delta-lognormal (two independent lognormals in-precip.).
  // Calculate precipitation fraction.
  if (l_use_precip_frac) {
    precip_fraction(workspace, gr, nzt, ngrdcol, hydromet_dim,                      // In
                    hydromet, cloud_frac, pdf_params_cloud_frac_1,                  // In
                    hm_metadata_l_mix_rat_hm, hm_metadata_l_frozen_hm,              // In
                    hm_metadata_hydromet_tol,                                       // In
                    pdf_params_cloud_frac_2, ice_supersat_frac,                     // In
                    pdf_params_ice_supersat_frac_1, pdf_params_ice_supersat_frac_2, // In
                    pdf_params_mixt_frac, clubb_params,                             // In
                    err_info,                                                       // InOut
                    precip_fracs_precip_frac,                                       // Out
                    precip_fracs_precip_frac_1,                                     // Out
                    precip_fracs_precip_frac_2,                                     // Out
                    precip_frac_tol,                                                // Out
                    stats                                                           // InOut
    );
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << " in setup_pdf_parameters_api after calling precip_fraction" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        return;
      }
    }
  } else {
    {
      const int slice_extent_0 = precip_fracs_precip_frac.extent_int(0);
      const int slice_extent_1 = precip_fracs_precip_frac.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_3", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                          = ((cell % slice_extent_0));
            const int slice_1                          = ((cell / slice_extent_0));
            precip_fracs_precip_frac(slice_0, slice_1) = one;
          });
    }
    {
      const int slice_extent_0 = precip_fracs_precip_frac_1.extent_int(0);
      const int slice_extent_1 = precip_fracs_precip_frac_1.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_4", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                            = ((cell % slice_extent_0));
            const int slice_1                            = ((cell / slice_extent_0));
            precip_fracs_precip_frac_1(slice_0, slice_1) = one;
          });
    }
    {
      const int slice_extent_0 = precip_fracs_precip_frac_2.extent_int(0);
      const int slice_extent_1 = precip_fracs_precip_frac_2.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_5", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                            = ((cell % slice_extent_0));
            const int slice_1                            = ((cell / slice_extent_0));
            precip_fracs_precip_frac_2(slice_0, slice_1) = one;
          });
    }
    {
      const int slice_extent_0 = precip_frac_tol.extent_int(0);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_6", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0        = cell;
            precip_frac_tol(slice_0) = cloud_frac_min;
          });
    }
  }
  // Calculate <N_cn> from Nc_in_cloud, whether Nc_in_cloud is predicted or
  // based on a prescribed value, and whether the value is constant or varying
  // over the grid level.
  if (!l_const_Nc_in_cloud) {
    // Ncn varies at each vertical level.
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::cells_7", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i                = ((cell % extent_1));
            const int k                = ((cell / extent_1));
            const_Ncnp2_on_Ncnm2(i, k) = hm_metadata_Ncnp2_on_Ncnm2;
          });
    }
    {
      const int slice_extent_0 = stdev_const_Ncnp2_on_Ncnm2.extent_int(0);
      const int slice_extent_1 = stdev_const_Ncnp2_on_Ncnm2.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_8", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            stdev_const_Ncnp2_on_Ncnm2(slice_0, slice_1) =
                stdev_L2N(const_Ncnp2_on_Ncnm2(slice_0, slice_1));
          });
    }
  } else {
    // l_const_Nc_in_cloud
    // Ncn is constant at each vertical level.
    {
      const int slice_extent_0 = const_Ncnp2_on_Ncnm2.extent_int(0);
      const int slice_extent_1 = const_Ncnp2_on_Ncnm2.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_9", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                      = ((cell % slice_extent_0));
            const int slice_1                      = ((cell / slice_extent_0));
            const_Ncnp2_on_Ncnm2(slice_0, slice_1) = zero;
          });
    }
    {
      const int slice_extent_0 = stdev_const_Ncnp2_on_Ncnm2.extent_int(0);
      const int slice_extent_1 = stdev_const_Ncnp2_on_Ncnm2.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_10", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                            = ((cell % slice_extent_0));
            const int slice_1                            = ((cell / slice_extent_0));
            stdev_const_Ncnp2_on_Ncnm2(slice_0, slice_1) = zero;
          });
    }
  }
  {
    const int slice_extent_0 = const_corr_chi_Ncn_n_cloud.extent_int(0);
    const int slice_extent_1 = const_corr_chi_Ncn_n_cloud.extent_int(1);
    Kokkos::parallel_for(
        "setup_pdf_parameters_api::array_11", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          const_corr_chi_Ncn_n_cloud(slice_0, slice_1) =
              corr_array_n_cloud(hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_chi);
        });
  }
  corr_NN2NL(nzt, ngrdcol,                                     // In
             const_corr_chi_Ncn_n_cloud,                       // In
             stdev_const_Ncnp2_on_Ncnm2, const_Ncnp2_on_Ncnm2, // In
             const_corr_chi_Ncn                                // Out
  );
  Nc_in_cloud_to_Ncnm(pdf_params_chi_1, pdf_params_chi_2, pdf_params_stdev_chi_1,
                      pdf_params_stdev_chi_2, pdf_params_mixt_frac, Nc_in_cloud,
                      pdf_params_cloud_frac_1, pdf_params_cloud_frac_2, const_Ncnp2_on_Ncnm2,
                      const_corr_chi_Ncn, Ncnm);
  // Calculate the overall variance of a precipitating hydrometeor (hm),
  // <hm'^2>.
  for (j = 0; j < hydromet_dim; ++j) {
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::cells_12", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            // Native loop: k = 1, nzt, 1
            if ((hydromet(i, k, j) >= hm_metadata_hydromet_tol(j))) {
              // There is some of the hydrometeor species found at level k.
              // Calculate the variance (overall) of the hydrometeor.
              hydrometp2_zt(i, k, j) =
                  ((((hm_metadata_hmp2_ip_on_hmm2_ip(j) + one) / precip_fracs_precip_frac(i, k)) -
                    one) *
                   (hydromet(i, k, j) * hydromet(i, k, j)));
            } else {
              hydrometp2_zt(i, k, j) = zero;
            }
          });
    }
  }
  // Interpolate the overall variance of a hydrometeor, <hm'^2>, to its home on
  // momentum grid levels.
  for (i = 0; i < hydromet_dim; ++i) {
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(hydrometp2_zt, Kokkos::ALL, Kokkos::ALL, i),
              Kokkos::subview(hydrometp2, Kokkos::ALL, Kokkos::ALL, i), zero_threshold);
    {
      const int slice_extent_0 = hydrometp2.extent_int(0);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_13", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                  = cell;
            hydrometp2(slice_0, gr_k_ub_zm, i) = zero;
          });
    }
  }
  // Calculate correlations involving w and Ncn by first calculating the
  // overall covariance of w and Ncn using the down-gradient approximation.
  if (l_calc_w_corr) {
    // Recalculate wm_zt and wp2_zt.  Mean vertical velocity may not be easy to
    // pass into this subroutine from a host model, and wp2_zt needs to have a
    // value consistent with the value it had when the PDF parameters involving w
    // were originally set in subroutine pdf_closure.  The variable wp2 has since
    // been advanced, resulting a new wp2_zt.  However, the value of wp2 here
    // needs to be consistent with wp2 at the time the PDF parameters were
    // calculated.
    // Calculate the overall mean of vertical velocity, w, on thermodynamic
    // levels.
    {
      const int slice_extent_0 = wm_zt.extent_int(0);
      const int slice_extent_1 = wm_zt.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_14", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0       = ((cell % slice_extent_0));
            const int slice_1       = ((cell / slice_extent_0));
            wm_zt(slice_0, slice_1) = compute_mean_binormal(pdf_params_w_1(slice_0, slice_1),
                                                            pdf_params_w_2(slice_0, slice_1),
                                                            pdf_params_mixt_frac(slice_0, slice_1));
          });
    }
    // Calculate the overall variance of vertical velocity on thermodynamic
    // levels.
    {
      const int slice_extent_0 = wp2_zt.extent_int(0);
      const int slice_extent_1 = wp2_zt.extent_int(1);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_15", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0        = ((cell % slice_extent_0));
            const int slice_1        = ((cell / slice_extent_0));
            wp2_zt(slice_0, slice_1) = compute_variance_binormal(
                wm_zt(slice_0, slice_1), pdf_params_w_1(slice_0, slice_1),
                pdf_params_w_2(slice_0, slice_1), sigma_w_1(slice_0, slice_1),
                sigma_w_2(slice_0, slice_1), pdf_params_mixt_frac(slice_0, slice_1));
          });
    }
    // Interpolate the covariances (overall) of w and precipitating
    // hydrometeors to thermodynamic grid levels.
    for (i = 0; i < hydromet_dim; ++i) {
      zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wphydrometp, Kokkos::ALL, Kokkos::ALL, i),
                Kokkos::subview(wphydrometp_zt, Kokkos::ALL, Kokkos::ALL, i));
    }
    for (j = 0; j < hydromet_dim; ++j) {
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "setup_pdf_parameters_api::cells_16", ((extent_0 * extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              if ((hydromet(i, k, j) < hm_metadata_hydromet_tol(j))) {
                wphydrometp_zt(i, k, j) = zero;
              }
            });
      }
    }
    // When the mean value of a precipitating hydrometeor is below tolerance
    // value, it is considered to have a value of 0, and the precipitating
    // hydrometeor does not vary over the grid level.  Any covariances
    // involving that precipitating hydrometeor also have values of 0 at that
    // grid level.
    for (j = 0; j < hydromet_dim; ++j) {
      // Clip the value of covariance <w'hm'> on thermodynamic levels.
      clip_covar(nzt, ngrdcol, clip_wphydrometp, wp2_zt,                         // In
                 Kokkos::subview(hydrometp2_zt, Kokkos::ALL, Kokkos::ALL, j),    // In
                 Kokkos::subview(wphydrometp_zt, Kokkos::ALL, Kokkos::ALL, j),   // InOut
                 Kokkos::subview(wphydrometp_chnge, Kokkos::ALL, Kokkos::ALL, j) // Out
      );
    }
    // Native loop: i = 1, hydromet_dim, 1
    {
      const int extent_0 = nzm;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::cells_17", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i        = ((cell % extent_1));
            const int k        = ((cell / extent_1));
            Kh_zm_c_K_hm(i, k) = (-(clubb_params(i, ic_K_hm)) * Kh_zm(i, k));
          });
    }
    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Kh_zm_c_K_hm, Ncnm,    // In
              wpNcnp_zm              // Out
    );
    // Boundary conditions; We are assuming zero flux.
    {
      const int slice_extent_0 = wpNcnp_zm.extent_int(0);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_18", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0     = cell;
            wpNcnp_zm(slice_0, 0) = zero;
          });
    }
    {
      const int slice_extent_0 = wpNcnp_zm.extent_int(0);
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::array_19", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = cell;
            wpNcnp_zm(slice_0, ((nzm - 1))) = zero;
          });
    }
    // Interpolate the covariances to thermodynamic grid levels.
    zm2zt_api(nzm, nzt, ngrdcol, gr, wpNcnp_zm, wpNcnp_zt);
    // When the mean value of Ncn is below tolerance value, it is considered
    // to have a value of 0, and Ncn does not vary over the grid level.  Any
    // covariance involving Ncn also has a value of 0 at that grid level.
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::cells_20", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            // Native loop: k = 1, nzt, 1
            if ((Ncnm(i, k) <= Ncn_tol)) {
              wpNcnp_zt(i, k) = zero;
            }
          });
    }
  }
  // l_calc_w_corr
  // Unpack CLUBB parameters
  {
    const int slice_extent_0 = omicron.extent_int(0);
    Kokkos::parallel_for(
        "setup_pdf_parameters_api::array_21", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = cell;
          omicron(slice_0)  = clubb_params(slice_0, iomicron);
        });
  }
  {
    const int slice_extent_0 = zeta_vrnce_rat.extent_int(0);
    Kokkos::parallel_for(
        "setup_pdf_parameters_api::array_22", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
          const int slice_0       = cell;
          zeta_vrnce_rat(slice_0) = clubb_params(slice_0, izeta_vrnce_rat);
        });
  }
  // !! Calculate the means and standard deviations involving PDF variables
  // !! -- w, chi, eta, N_cn, and any precipitating hydrometeors (hm in-precip)
  // !! -- for each PDF component.
  compute_mean_stdev(nzt, ngrdcol, hydromet_dim, // In
                     hydromet, hydrometp2_zt,    // In
                     hm_metadata,                // In
                     Ncnm, pdf_params_mixt_frac, // In
                     precip_fracs_precip_frac,   // In
                     precip_fracs_precip_frac_1, // In
                     precip_fracs_precip_frac_2, // In
                     precip_frac_tol,            // In
                     pdf_params_w_1,             // In
                     pdf_params_w_2,             // In
                     sigma_w_1,                  // In
                     sigma_w_1,                  // In
                     pdf_params_chi_1,           // In
                     pdf_params_chi_2,           // In
                     pdf_params_stdev_chi_1,     // In
                     pdf_params_stdev_chi_2,     // In
                     pdf_params_stdev_eta_1,     // In
                     pdf_params_stdev_eta_2,     // In
                     pdf_params_thl_1,           // In
                     pdf_params_thl_2,           // In
                     pdf_dim,                    // In
                     omicron, zeta_vrnce_rat,    // In
                     l_const_Nc_in_cloud,        // In
                     mu_x_1, mu_x_2,             // Out
                     sigma_x_1, sigma_x_2,       // Out
                     hm_1, hm_2,                 // Out
                     sigma2_on_mu2_ip_1,         // Out
                     sigma2_on_mu2_ip_2          // Out
  );
  // !! Transform the component means and standard deviations involving
  // !! precipitating hydrometeors (hm in-precip) and N_cn -- ln hm and
  // !! ln N_cn -- to normal space for each PDF component.
  norm_transform_mean_stdev(nzt, ngrdcol,              // In
                            hydromet_dim, hm_metadata, // In
                            hm_1, hm_2,                // In
                            Ncnm, pdf_dim,             // In
                            mu_x_1, mu_x_2,            // In
                            sigma_x_1, sigma_x_2,      // In
                            sigma2_on_mu2_ip_1,        // In
                            sigma2_on_mu2_ip_2,        // In
                            l_const_Nc_in_cloud,       // In
                            mu_x_1_n, mu_x_2_n,        // Out
                            sigma_x_1_n, sigma_x_2_n   // Out
  );
  // !! Calculate the normal space correlations.
  // !! The normal space correlations are the same as the true correlations
  // !! except when at least one of the variables involved is a precipitating
  // !! hydrometeor or Ncn.  In these cases, the normal space correlation
  // !! involves the natural logarithm of the precipitating hydrometeors, ln hm
  // !! (for example, ln r_r and ln N_r), and ln N_cn for each PDF component.
  // The lower supported diagnose API rejects l_calc_w_corr=true: Native leaves
  // its working matrix undefined in that conjunction, as recorded by JAX.
  if (l_diagnose_correlations) {
    const auto rcm_pdf_host = to_host_space(workspace, rcm_pdf);
    for (k = 0; k < nzt; ++k) {
      for (i = 0; i < ngrdcol; ++i) {
        if ((rcm_pdf_host(i, k) > rc_tol)) {
          diagnose_correlations(
              workspace, pdf_dim, hm_metadata_iiPDF_w,                        // In
              corr_array_n_cloud,                                             // In
              l_calc_w_corr,                                                  // In
              Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL) // Out
          );
          diagnose_correlations(
              workspace, pdf_dim, hm_metadata_iiPDF_w,                        // In
              corr_array_n_cloud,                                             // In
              l_calc_w_corr,                                                  // In
              Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL) // Out
          );
        } else {
          diagnose_correlations(
              workspace, pdf_dim, hm_metadata_iiPDF_w,                        // In
              corr_array_n_below,                                             // In
              l_calc_w_corr,                                                  // In
              Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL) // Out
          );
          diagnose_correlations(
              workspace, pdf_dim, hm_metadata_iiPDF_w,                        // In
              corr_array_n_below,                                             // In
              l_calc_w_corr,                                                  // In
              Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL) // Out
          );
        }
      }
    }
    for (k = 0; k < nzt; ++k) {
      for (i = 0; i < ngrdcol; ++i) {
        calc_cholesky_corr_mtx_approx(
            workspace, pdf_dim, hm_metadata_iiPDF_w,                              // In
            Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL),      // In
            Kokkos::subview(corr_cholesky_mtx_1, i, k, Kokkos::ALL, Kokkos::ALL), // Out
            Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL)       // Out
        );
      }
    }
    for (k = 0; k < nzt; ++k) {
      for (i = 0; i < ngrdcol; ++i) {
        calc_cholesky_corr_mtx_approx(
            workspace, pdf_dim, hm_metadata_iiPDF_w,                              // In
            Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL),      // In
            Kokkos::subview(corr_cholesky_mtx_2, i, k, Kokkos::ALL, Kokkos::ALL), // Out
            Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL)       // Out
        );
      }
    }
  } else {
    // if .not. l_diagnose_correlations
    if ((l_fix_w_chi_eta_correlations && !l_calc_w_corr)) {
      // When the flags are set this way, the correlation matrices do not vary with any vertical
      // values, and instead are determined entirely by prescribed values. This results in there
      // being only two unique correlation matrices, one for when the grid box is in cloud
      // and one for when it is not. So instead of setting up correlation matrices for all
      // grid boxes then calculating their Cholesky decomps, we can simply set up two correlation
      // matrices, one for in cloud and one for out cloud, calculate the corresponding
      // Cholesky decompositions, then use the value of rc at each grid box to determine whether
      // we assign the in cloud or out of cloud matrices to that grid box.
      calc_corr_norm_and_cholesky_factor(workspace, nzt, ngrdcol, iiPDF_type,     // In
                                         pdf_dim, hm_metadata,                    // In
                                         pdf_params_rc_1, pdf_params_rc_2,        // In
                                         corr_array_n_cloud, corr_array_n_below,  // In
                                         corr_array_1_n, corr_array_2_n,          // Out
                                         corr_cholesky_mtx_1, corr_cholesky_mtx_2 // Out
      );
    } else {
      // The correlation matrices can vary with vertical values. So we need to set the
      // correlation matrices up for each grid box, then find the Cholesky decomp for each
      // grid box individually. This is very computationally expensive.
      comp_corr_norm(workspace, nzt, pdf_dim, ngrdcol, hydromet_dim, hm_metadata, // In
                     wm_zt, pdf_params_rc_1, pdf_params_rc_2,                     // In
                     pdf_params_mixt_frac,                                        // In
                     precip_fracs_precip_frac_1,                                  // In
                     precip_fracs_precip_frac_2,                                  // In
                     wpNcnp_zt, wphydrometp_zt,                                   // In
                     mu_x_1, mu_x_2,                                              // In
                     sigma_x_1, sigma_x_2,                                        // In
                     sigma_x_1_n, sigma_x_2_n,                                    // In
                     corr_array_n_cloud, corr_array_n_below,                      // In
                     pdf_params,                                                  // In
                     iiPDF_type,                                                  // In
                     l_calc_w_corr,                                               // In
                     l_fix_w_chi_eta_correlations,                                // In
                     corr_array_1_n, corr_array_2_n                               // Out
      );
      // Compute choleksy factorization for the correlation matrix of 1st PDF component
      for (k = 0; k < nzt; ++k) {
        for (i = 0; i < ngrdcol; ++i) {
          Cholesky_factor(
              workspace,
              pdf_dim,                                                              // In
              Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL),      // In
              corr_array_scaling,                                                   // Out
              Kokkos::subview(corr_cholesky_mtx_1, i, k, Kokkos::ALL, Kokkos::ALL), // Out
              l_corr_array_scaling                                                  // Out
          );
        }
      }
      // Compute choleksy factorization for the correlation matrix of 2nd PDF component
      for (k = 0; k < nzt; ++k) {
        for (i = 0; i < ngrdcol; ++i) {
          Cholesky_factor(
              workspace,
              pdf_dim,                                                              // In
              Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL),      // In
              corr_array_scaling,                                                   // Out
              Kokkos::subview(corr_cholesky_mtx_2, i, k, Kokkos::ALL, Kokkos::ALL), // Out
              l_corr_array_scaling                                                  // Out
          );
        }
      }
    }
  }
  // l_diagnose_correlations
  // hydromet_pdf_params is optional, so if it is not present we simply skip the steps
  // to compute it.
  if (hydromet_pdf_params != nullptr) {
    // !! Calculate the true correlations for each PDF component.
    denorm_transform_corr(nzt, ngrdcol, pdf_dim, hm_metadata,     // In
                          sigma_x_1_n, sigma_x_2_n,               // In
                          sigma2_on_mu2_ip_1, sigma2_on_mu2_ip_2, // In
                          corr_array_1_n,                         // In
                          corr_array_2_n,                         // In
                          corr_array_1, corr_array_2              // Out
    );
    // !! Pack the PDF parameters
    pack_hydromet_pdf_params(nzt, ngrdcol,               // In
                             hydromet_dim, hm_metadata,  // In
                             hm_1, hm_2, pdf_dim,        // In
                             mu_x_1, mu_x_2,             // In
                             sigma_x_1, sigma_x_2,       // In
                             corr_array_1, corr_array_2, // In
                             *hydromet_pdf_params        // Out
    );
  }
  if (stats.l_sample) {
    for (j = 0; j < hydromet_dim; ++j) {
      hm_type  = std::string(hm_metadata_hydromet_list_host.data() +
                                 ((hm_metadata_hydromet_list_host.extent(0) * j)),
                             hm_metadata_hydromet_list_host.extent(0));
      var_name = setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "p2_zt";
      stats_update(var_name,                                                    // In
                   Kokkos::subview(hydrometp2_zt, Kokkos::ALL, Kokkos::ALL, j), // In
                   stats                                                        // InOut
      );
    }
    // TODO: Native samples true-space correlation scratch even when the optional
    // hydrometeor output is absent. That conjunction is undefined in Native;
    // retain the call/guard order without adopting the different JAX guard.
    pdf_param_hm_stats(nzt, ngrdcol, pdf_dim, hydromet_dim, hm_metadata, // In
                       hm_1, hm_2,                                       // In
                       mu_x_1, mu_x_2,                                   // In
                       sigma_x_1, sigma_x_2,                             // In
                       corr_array_1, corr_array_2,                       // In
                       stats                                             // InOut
    );
    // !! Statistics for normal space PDF parameters involving hydrometeors.
    pdf_param_ln_hm_stats(workspace, nzt, ngrdcol, pdf_dim, hm_metadata, // In
                          mu_x_1_n, mu_x_2_n,                            // In
                          sigma_x_1_n, sigma_x_2_n,                      // In
                          corr_array_1_n, corr_array_2_n,                // In
                          stats                                          // InOut
    );
    if (var_on_stats_list(stats, "rtp2_from_chi")) {
      // Preserve Native's post-DO j=hydromet_dim+1 actual rt_2 column.
      // TODO: Native reads outside rt_2 when hydromet_dim+1>ngrdcol. No corrected
      // JAX indexing or fabricated value is substituted at this reference boundary.
      const int rt_2_column = hydromet_dim;
      Kokkos::parallel_for(
          "setup_pdf_parameters_api::rtp2_from_chi", ((ngrdcol * nzt)),
          KOKKOS_LAMBDA(const int cell) {
            const int i            = ((cell % ngrdcol));
            const int k            = ((cell / ngrdcol));
            rtp2_zt_from_chi(i, k) = compute_rtp2_from_chi(
                pdf_params_stdev_chi_1(i, k), pdf_params_stdev_chi_2(i, k),
                pdf_params_stdev_eta_1(i, k), pdf_params_stdev_eta_2(i, k), pdf_params_rt_1(i, k),
                pdf_params_rt_2(rt_2_column, k), pdf_params_crt_1(i, k), pdf_params_crt_2(i, k),
                pdf_params_mixt_frac(i, k),
                corr_array_1_n(i, k, hm_metadata_iiPDF_chi, hm_metadata_iiPDF_eta),
                corr_array_2_n(i, k, hm_metadata_iiPDF_chi, hm_metadata_iiPDF_eta));
          });
      zt2zm_api(nzm, nzt, ngrdcol, gr, rtp2_zt_from_chi, rtp2_zm_from_chi, zero_threshold);
      stats_update("rtp2_from_chi", rtp2_zm_from_chi, // In
                   stats                              // InOut
      );
    }
    stats_update("precip_frac", precip_fracs_precip_frac, // In
                 stats                                    // InOut
    );
    stats_update("precip_frac_1", precip_fracs_precip_frac_1, // In
                 stats                                        // InOut
    );
    stats_update("precip_frac_2", precip_fracs_precip_frac_2, // In
                 stats                                        // InOut
    );
    stats_update("Ncnm", Ncnm, // In
                 stats         // InOut
    );
  }
  if (clubb_at_least_debug_level_api(2)) {
    for (k = 0; k < nzt; ++k) {
      for (i = 0; i < ngrdcol; ++i) {
        assert_corr_symmetric(workspace,
                              Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL), // In
                              pdf_dim,                                                         // In
                              err_info // InOut
        );
        assert_corr_symmetric(workspace,
                              Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL), // In
                              pdf_dim,                                                         // In
                              err_info // InOut
        );
        if ((to_host_space(workspace, err_info.err_code)(i) == clubb_fatal_error)) {
          {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            const auto header_host = to_host_space(workspace, err_info.err_header);
            diagnostic.write(header_host.data() + ((header_host.extent(0) * i)),
                             header_host.extent(0));
            std::cerr << diagnostic.str() << '\n';
          }
        }
      }
    }
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << " in setup_pdf_parameters_api after calling assert_corr_symmetric" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
    }
  }
  return;
}
// =============================================================================
void compute_mean_stdev(int nzt, int ngrdcol, int hydromet_dim, // In
                        Real3D hydromet, Real3D hydrometp2_zt,  // In
                        const hm_metadata_type& hm_metadata,    // In
                        Real2D Ncnm, Real2D mixt_frac,          // In
                        Real2D precip_frac,                     // In
                        Real2D precip_frac_1,                   // In
                        Real2D precip_frac_2,                   // In
                        Real1D precip_frac_tol,                 // In
                        Real2D w_1, Real2D w_2,                 // In
                        Real2D stdev_w_1, Real2D stdev_w_2,     // In
                        Real2D chi_1, Real2D chi_2,             // In
                        Real2D stdev_chi_1, Real2D stdev_chi_2, // In
                        Real2D stdev_eta_1, Real2D stdev_eta_2, // In
                        Real2D thl_1, Real2D thl_2,             // In
                        int pdf_dim,                            // In
                        Real1D omicron, Real1D zeta_vrnce_rat,  // In
                        bool l_const_Nc_in_cloud,               // In
                        Real3D mu_x_1, Real3D mu_x_2,           // Out
                        Real3D sigma_x_1, Real3D sigma_x_2,     // Out
                        Real3D hm_1, Real3D hm_2,               // Out
                        Real3D sigma_hm_1_sqd_on_mu_hm_1_sqd,   // Out
                        Real3D sigma_hm_2_sqd_on_mu_hm_2_sqd    // Out
) {
  // Description:
  // Calculates the means and standard deviations (for each PDF component) of
  // chi, eta, w, Ncn, and the precipitating hydrometeors.  For the
  // precipitating hydrometeors, the component means and standard deviations
  // are in-precip.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // hydromet: Mean of a hydrometeor (overall)         [hm units]
  // hydrometp2_zt: Variance of a hydrometeor (overall)     [(hm units)^2]
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // Ncnm: Mean simplified cloud nuclei concentration   [num/kg]
  // mixt_frac: Mixture fraction                             [-]
  // precip_frac: Precipitation fraction (overall)             [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)   [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)   [-]
  // precip_frac_tol: Minimum precip. frac. when hydromet. are present  [-]
  // w_1: Mean of w (1st PDF component)                  [m/s]
  // w_2: Mean of w (2nd PDF component)                  [m/s]
  // stdev_w_1: Standard deviation of w (1st PDF component)              [m/s]
  // stdev_w_2: Standard deviation of w (2nd PDF component)              [m/s]
  // chi_1: Mean of chi (1st PDF component)                [kg/kg]
  // chi_2: Mean of chi (2nd PDF component)                [kg/kg]
  // stdev_chi_1: Standard deviation of chi (1st PDF component)  [kg/kg]
  // stdev_chi_2: Standard deviation of chi (2nd PDF component)  [kg/kg]
  // stdev_eta_1: Standard deviation of eta (1st PDF component)  [kg/kg]
  // stdev_eta_2: Standard deviation of eta (2nd PDF component)  [kg/kg]
  // thl_1: Mean of thl (1st PDF component)                [K]
  // thl_2: Mean of thl (2nd PDF component)                [K]
  // pdf_dim: Number of PDF variables
  // omicron: Relative width parameter, omicron = R / Rmax    [-]
  // zeta_vrnce_rat: Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]
  // l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K)
  // mu_x_1: Mean array of PDF vars. (1st PDF component)    [units vary]
  // mu_x_2: Mean array of PDF vars. (2nd PDF component)    [units vary]
  // sigma_x_1: Standard deviation array of PDF vars (comp. 1) [units vary]
  // sigma_x_2: Standard deviation array of PDF vars (comp. 2) [units vary]
  // hm_1: Mean of a precip. hydrometeor (1st PDF component)   [units vary]
  // hm_2: Mean of a precip. hydrometeor (2nd PDF component)   [units vary]
  // sigma_hm_1_sqd_on_mu_hm_1_sqd: Ratio sigma_hm_1^2 / mu_hm_1^2      [-]
  // sigma_hm_2_sqd_on_mu_hm_2_sqd: Ratio sigma_hm_2^2 / mu_hm_2^2      [-]

  int ivar{};
  int hm_idx{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_w                 = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_chi               = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta               = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_Ncn               = hm_metadata.iiPDF_Ncn;
  const auto hm_metadata_Ncnp2_on_Ncnm2          = hm_metadata.Ncnp2_on_Ncnm2;
  const auto hm_metadata_hmp2_ip_on_hmm2_ip      = hm_metadata.hmp2_ip_on_hmm2_ip;
  const auto hm_metadata_hydromet_tol            = hm_metadata.hydromet_tol;
  const auto hm_metadata_hmp2_ip_on_hmm2_ip_host = to_host_space(hm_metadata_hmp2_ip_on_hmm2_ip);
  const auto hm_metadata_hydromet_tol_host       = to_host_space(hm_metadata_hydromet_tol);

  // !! Initialize output variables.
  {
    const int slice_extent_0 = hm_1.extent_int(0);
    const int slice_extent_1 = hm_1.extent_int(1);
    const int slice_extent_2 = hm_1.extent_int(2);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_1", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0               = ((cell % slice_extent_0));
          const int slice_1               = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2               = ((cell / (slice_extent_0 * slice_extent_1)));
          hm_1(slice_0, slice_1, slice_2) = zero;
        });
  }
  {
    const int slice_extent_0 = hm_2.extent_int(0);
    const int slice_extent_1 = hm_2.extent_int(1);
    const int slice_extent_2 = hm_2.extent_int(2);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_2", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0               = ((cell % slice_extent_0));
          const int slice_1               = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2               = ((cell / (slice_extent_0 * slice_extent_1)));
          hm_2(slice_0, slice_1, slice_2) = zero;
        });
  }
  {
    const int slice_extent_0 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(0);
    const int slice_extent_1 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(1);
    const int slice_extent_2 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(2);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_3", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2 = ((cell / (slice_extent_0 * slice_extent_1)));
          sigma_hm_1_sqd_on_mu_hm_1_sqd(slice_0, slice_1, slice_2) = zero;
        });
  }
  {
    const int slice_extent_0 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(0);
    const int slice_extent_1 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(1);
    const int slice_extent_2 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(2);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_4", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2 = ((cell / (slice_extent_0 * slice_extent_1)));
          sigma_hm_2_sqd_on_mu_hm_2_sqd(slice_0, slice_1, slice_2) = zero;
        });
  }
  // !! Enter the PDF parameters.
  // !! Vertical velocity, w.
  // Mean of vertical velocity, w, in PDF component 1.
  {
    const int slice_extent_0 = mu_x_1.extent_int(0);
    const int slice_extent_1 = mu_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_5", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                             = ((cell % slice_extent_0));
          const int slice_1                             = ((cell / slice_extent_0));
          mu_x_1(slice_0, slice_1, hm_metadata_iiPDF_w) = w_1(slice_0, slice_1);
        });
  }
  // Mean of vertical velocity, w, in PDF component 2.
  {
    const int slice_extent_0 = mu_x_2.extent_int(0);
    const int slice_extent_1 = mu_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_6", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                             = ((cell % slice_extent_0));
          const int slice_1                             = ((cell / slice_extent_0));
          mu_x_2(slice_0, slice_1, hm_metadata_iiPDF_w) = w_2(slice_0, slice_1);
        });
  }
  // Standard deviation of vertical velocity, w, in PDF component 1.
  {
    const int slice_extent_0 = sigma_x_1.extent_int(0);
    const int slice_extent_1 = sigma_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_7", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                = ((cell % slice_extent_0));
          const int slice_1                                = ((cell / slice_extent_0));
          sigma_x_1(slice_0, slice_1, hm_metadata_iiPDF_w) = stdev_w_1(slice_0, slice_1);
        });
  }
  // Standard deviation of vertical velocity, w, in PDF component 2.
  {
    const int slice_extent_0 = sigma_x_2.extent_int(0);
    const int slice_extent_1 = sigma_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_8", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                = ((cell % slice_extent_0));
          const int slice_1                                = ((cell / slice_extent_0));
          sigma_x_2(slice_0, slice_1, hm_metadata_iiPDF_w) = stdev_w_2(slice_0, slice_1);
        });
  }
  // !! Extended liquid water mixing ratio, chi.
  // Mean of extended liquid water mixing ratio, chi (old s),
  // in PDF component 1.
  {
    const int slice_extent_0 = mu_x_1.extent_int(0);
    const int slice_extent_1 = mu_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_9", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_1(slice_0, slice_1, hm_metadata_iiPDF_chi) = chi_1(slice_0, slice_1);
        });
  }
  // Mean of extended liquid water mixing ratio, chi (old s),
  // in PDF component 2.
  {
    const int slice_extent_0 = mu_x_2.extent_int(0);
    const int slice_extent_1 = mu_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_10", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_2(slice_0, slice_1, hm_metadata_iiPDF_chi) = chi_2(slice_0, slice_1);
        });
  }
  // Standard deviation of extended liquid water mixing ratio, chi (old s),
  // in PDF component 1.
  {
    const int slice_extent_0 = sigma_x_1.extent_int(0);
    const int slice_extent_1 = sigma_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_11", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                  = ((cell % slice_extent_0));
          const int slice_1                                  = ((cell / slice_extent_0));
          sigma_x_1(slice_0, slice_1, hm_metadata_iiPDF_chi) = stdev_chi_1(slice_0, slice_1);
        });
  }
  // Standard deviation of extended liquid water mixing ratio, chi (old s),
  // in PDF component 2.
  {
    const int slice_extent_0 = sigma_x_2.extent_int(0);
    const int slice_extent_1 = sigma_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_12", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                  = ((cell % slice_extent_0));
          const int slice_1                                  = ((cell / slice_extent_0));
          sigma_x_2(slice_0, slice_1, hm_metadata_iiPDF_chi) = stdev_chi_2(slice_0, slice_1);
        });
  }
  // !! Coordinate orthogonal to chi, eta.
  // Mean of eta (old t) in PDF component 1.
  // Set the component mean values of eta to 0.
  // The component mean values of eta are not important.  They can be set to
  // anything.  They cancel out in the model code.  However, the best thing to
  // do is to set them to 0 and avoid any kind of numerical error.
  {
    const int slice_extent_0 = mu_x_1.extent_int(0);
    const int slice_extent_1 = mu_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_13", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_1(slice_0, slice_1, hm_metadata_iiPDF_eta) = zero;
        });
  }
  // Mean of eta (old t) in PDF component 2.
  // Set the component mean values of eta to 0.
  // The component mean values of eta are not important.  They can be set to
  // anything.  They cancel out in the model code.  However, the best thing to
  // do is to set them to 0 and avoid any kind of numerical error.
  {
    const int slice_extent_0 = mu_x_2.extent_int(0);
    const int slice_extent_1 = mu_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_14", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_2(slice_0, slice_1, hm_metadata_iiPDF_eta) = zero;
        });
  }
  // Standard deviation of eta (old t) in PDF component 1.
  {
    const int slice_extent_0 = sigma_x_1.extent_int(0);
    const int slice_extent_1 = sigma_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_15", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                  = ((cell % slice_extent_0));
          const int slice_1                                  = ((cell / slice_extent_0));
          sigma_x_1(slice_0, slice_1, hm_metadata_iiPDF_eta) = stdev_eta_1(slice_0, slice_1);
        });
  }
  // Standard deviation of eta (old t) in PDF component 2.
  {
    const int slice_extent_0 = sigma_x_2.extent_int(0);
    const int slice_extent_1 = sigma_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_16", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                                  = ((cell % slice_extent_0));
          const int slice_1                                  = ((cell / slice_extent_0));
          sigma_x_2(slice_0, slice_1, hm_metadata_iiPDF_eta) = stdev_eta_2(slice_0, slice_1);
        });
  }
  // !! Simplified cloud nuclei concentration, Ncn.
  // Mean of simplified cloud nuclei concentration, Ncn, in PDF component 1.
  {
    const int slice_extent_0 = mu_x_1.extent_int(0);
    const int slice_extent_1 = mu_x_1.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_17", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_1(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = Ncnm(slice_0, slice_1);
        });
  }
  // Mean of simplified cloud nuclei concentration, Ncn, in PDF component 2.
  {
    const int slice_extent_0 = mu_x_2.extent_int(0);
    const int slice_extent_1 = mu_x_2.extent_int(1);
    Kokkos::parallel_for(
        "compute_mean_stdev::array_18", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                               = ((cell % slice_extent_0));
          const int slice_1                               = ((cell / slice_extent_0));
          mu_x_2(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = Ncnm(slice_0, slice_1);
        });
  }
  // Standard deviation of simplified cloud nuclei concentration, Ncn,
  // in PDF component 1.
  if (!l_const_Nc_in_cloud) {
    // Ncn varies in both PDF components.
    {
      const int slice_extent_0 = sigma_x_1.extent_int(0);
      const int slice_extent_1 = sigma_x_1.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_19", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_1(slice_0, slice_1, hm_metadata_iiPDF_Ncn) =
                (Kokkos::sqrt(hm_metadata_Ncnp2_on_Ncnm2) * Ncnm(slice_0, slice_1));
          });
    }
    {
      const int slice_extent_0 = sigma_x_2.extent_int(0);
      const int slice_extent_1 = sigma_x_2.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_20", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_2(slice_0, slice_1, hm_metadata_iiPDF_Ncn) =
                (Kokkos::sqrt(hm_metadata_Ncnp2_on_Ncnm2) * Ncnm(slice_0, slice_1));
          });
    }
    // Ncn is not an official hydrometeor.  However, both the
    // sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
    // have size pdf_dim, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
    // sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
    {
      const int slice_extent_0 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(0);
      const int slice_extent_1 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_21", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_hm_1_sqd_on_mu_hm_1_sqd(slice_0, slice_1, hm_metadata_iiPDF_Ncn) =
                hm_metadata_Ncnp2_on_Ncnm2;
          });
    }
    {
      const int slice_extent_0 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(0);
      const int slice_extent_1 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_22", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_hm_2_sqd_on_mu_hm_2_sqd(slice_0, slice_1, hm_metadata_iiPDF_Ncn) =
                hm_metadata_Ncnp2_on_Ncnm2;
          });
    }
  } else {
    // l_const_Nc_in_cloud
    // Ncn is constant in both PDF components.
    {
      const int slice_extent_0 = sigma_x_1.extent_int(0);
      const int slice_extent_1 = sigma_x_1.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_23", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                  = ((cell % slice_extent_0));
            const int slice_1                                  = ((cell / slice_extent_0));
            sigma_x_1(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = zero;
          });
    }
    {
      const int slice_extent_0 = sigma_x_2.extent_int(0);
      const int slice_extent_1 = sigma_x_2.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_24", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                  = ((cell % slice_extent_0));
            const int slice_1                                  = ((cell / slice_extent_0));
            sigma_x_2(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = zero;
          });
    }
    // Ncn is not an official hydrometeor.  However, both the
    // sigma_hm_1_sqd_on_mu_hm_1_sqd and sigma_hm_2_sqd_on_mu_hm_2_sqd arrays
    // have size pdf_dim, and both sigma_Ncn_1^2/mu_Ncn_1^2 and
    // sigma_Ncn_2^2/mu_Ncn_2^2 need to be output as part of these arrays.
    {
      const int slice_extent_0 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(0);
      const int slice_extent_1 = sigma_hm_1_sqd_on_mu_hm_1_sqd.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_25", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_hm_1_sqd_on_mu_hm_1_sqd(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = zero;
          });
    }
    {
      const int slice_extent_0 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(0);
      const int slice_extent_1 = sigma_hm_2_sqd_on_mu_hm_2_sqd.extent_int(1);
      Kokkos::parallel_for(
          "compute_mean_stdev::array_26", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_hm_2_sqd_on_mu_hm_2_sqd(slice_0, slice_1, hm_metadata_iiPDF_Ncn) = zero;
          });
    }
  }
  // .not. l_const_Nc_in_cloud
  // !! Precipitating hydrometeor species.
  for (ivar = ((hm_metadata_iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
    hm_idx = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
    calc_comp_mu_sigma_hm(
        nzt, ngrdcol,                                                                   // In
        Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, hm_idx),                    // In
        Kokkos::subview(hydrometp2_zt, Kokkos::ALL, Kokkos::ALL, hm_idx),               // In
        hm_metadata_hmp2_ip_on_hmm2_ip_host(hm_idx),                                    // In
        mixt_frac, precip_frac,                                                         // In
        precip_frac_1, precip_frac_2,                                                   // In
        hm_metadata_hydromet_tol_host(hm_idx),                                          // In
        precip_frac_tol,                                                                // In
        thl_1, thl_2,                                                                   // In
        omicron, zeta_vrnce_rat,                                                        // In
        Kokkos::subview(mu_x_1, Kokkos::ALL, Kokkos::ALL, ivar),                        // Out
        Kokkos::subview(mu_x_2, Kokkos::ALL, Kokkos::ALL, ivar),                        // Out
        Kokkos::subview(sigma_x_1, Kokkos::ALL, Kokkos::ALL, ivar),                     // Out
        Kokkos::subview(sigma_x_2, Kokkos::ALL, Kokkos::ALL, ivar),                     // Out
        Kokkos::subview(hm_1, Kokkos::ALL, Kokkos::ALL, hm_idx),                        // Out
        Kokkos::subview(hm_2, Kokkos::ALL, Kokkos::ALL, hm_idx),                        // Out
        Kokkos::subview(sigma_hm_1_sqd_on_mu_hm_1_sqd, Kokkos::ALL, Kokkos::ALL, ivar), // Out
        Kokkos::subview(sigma_hm_2_sqd_on_mu_hm_2_sqd, Kokkos::ALL, Kokkos::ALL, ivar)  // Out
    );
  }
  // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
  return;
}
// =============================================================================
namespace {
void calc_corr_norm_and_cholesky_factor(ScratchWorkspace& workspace,                      // InOut
                                        int nzt, int ngrdcol, int iiPDF_type,             // In
                                        int pdf_dim, const hm_metadata_type& hm_metadata, // In
                                        Real2D rc_1, Real2D rc_2,                         // In
                                        Real2D corr_array_n_cloud, Real2D corr_array_n_below, // In
                                        Real4D corr_array_1_n, Real4D corr_array_2_n,         // Out
                                        Real4D corr_cholesky_mtx_1,
                                        Real4D corr_cholesky_mtx_2 // Out
) {
  // Description: This subroutine computes the correlation arrays and correlation
  //   Cholesky matrices of PDF vars for both components. Here, we assume that
  //   there are only two unique correlation arrays, which allows us to compute
  //   these two unique arrays and their corresponding Cholesky decompositions,
  //   then use rc to determine which one to assign to each grid column and
  //   vertical level. If the correlation arrays vary based on vertically varying
  //   values, then this subroutine is not appropriate.
  //
  // References:
  //   https://github.com/larson-group/cam/issues/129#issuecomment-816205563
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical levels
  // ngrdcol: Number of grid columns
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF.
  // pdf_dim: Number of variables in the corr/mean/stdev arrays
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // rc_1: Mean of r_c (1st PDF component)                 [kg/kg]
  // rc_2: Mean of r_c (2nd PDF component)                 [kg/kg]
  // corr_array_n_cloud: Prescribed correlation array in cloud        [-]
  // corr_array_n_below: Prescribed correlation array below cloud     [-]
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // corr_cholesky_mtx_1: Transposed corr. cholesky matrix, 1st comp.        [-]
  // corr_cholesky_mtx_2: Transposed corr. cholesky matrix, 2nd comp.        [-]

  auto workspace_scope = workspace.scope();
  Real2D corr_array_cloud =
      workspace_scope.view<Real2D>("calc_corr_norm_and_cholesky_factor::corr_array_cloud", pdf_dim,
                                   pdf_dim); // General in cloud corr. matrix
  Real2D corr_array_below =
      workspace_scope.view<Real2D>("calc_corr_norm_and_cholesky_factor::corr_array_below", pdf_dim,
                                   pdf_dim); // General out of cloud corr. matrix
  Real2D corr_cholesky_mtx_cloud =
      workspace_scope.view<Real2D>("calc_corr_norm_and_cholesky_factor::corr_cholesky_mtx_cloud",
                                   pdf_dim, pdf_dim); // General in cloud Cholesky matrix
  Real2D corr_cholesky_mtx_below =
      workspace_scope.view<Real2D>("calc_corr_norm_and_cholesky_factor::corr_cholesky_mtx_below",
                                   pdf_dim, pdf_dim); // General out of cloud Cholesky matrix
  bool l_corr_array_scaling{}; // Dummy variable that we need for calling Cholesky_factor
  Real1D corr_array_scaling = workspace_scope.view<Real1D>(
      "calc_corr_norm_and_cholesky_factor::corr_array_scaling",
      pdf_dim); // Dummy variable that we need for calling Cholesky_factor
  int ivar{};
  int jvar{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_w   = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_chi = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_Ncn = hm_metadata.iiPDF_Ncn;

  // Initialize correlation arrays with prescribed values
  {
    const int slice_extent_0 = corr_array_cloud.extent_int(0);
    const int slice_extent_1 = corr_array_cloud.extent_int(1);
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::array_1", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                  = ((cell % slice_extent_0));
          const int slice_1                  = ((cell / slice_extent_0));
          corr_array_cloud(slice_0, slice_1) = corr_array_n_cloud(slice_0, slice_1);
        });
  }
  {
    const int slice_extent_0 = corr_array_below.extent_int(0);
    const int slice_extent_1 = corr_array_below.extent_int(1);
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::array_2", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                  = ((cell % slice_extent_0));
          const int slice_1                  = ((cell / slice_extent_0));
          corr_array_below(slice_0, slice_1) = corr_array_n_below(slice_0, slice_1);
        });
  }
  // The ADG1 PDF fixes the correlation of w and rt and the correlation of
  // w and theta_l to be 0, which means the correlation of w and chi and the
  // correlation of w and eta must also be 0.
  if ((((iiPDF_type == iiPDF_ADG1) || (iiPDF_type == iiPDF_ADG2) ||
        (iiPDF_type == iiPDF_new_hybrid)) &&
       l_follow_ADG1_PDF_standards)) {
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_3", 1, KOKKOS_LAMBDA(const int) {
          corr_array_cloud(hm_metadata_iiPDF_w, hm_metadata_iiPDF_chi) = zero;
        });
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_4", 1, KOKKOS_LAMBDA(const int) {
          corr_array_below(hm_metadata_iiPDF_w, hm_metadata_iiPDF_chi) = zero;
        });
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_5", 1, KOKKOS_LAMBDA(const int) {
          corr_array_cloud(hm_metadata_iiPDF_w, hm_metadata_iiPDF_eta) = zero;
        });
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_6", 1, KOKKOS_LAMBDA(const int) {
          corr_array_below(hm_metadata_iiPDF_w, hm_metadata_iiPDF_eta) = zero;
        });
  }
  // Ncn is an inherently in-cloud property, so replace out of cloud correlation values
  // with in cloud ones.
  Kokkos::parallel_for(
      "calc_corr_norm_and_cholesky_factor::scalar_7", 1, KOKKOS_LAMBDA(const int) {
        corr_array_below(hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_chi) =
            corr_array_cloud(hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_chi);
      });
  Kokkos::parallel_for(
      "calc_corr_norm_and_cholesky_factor::scalar_8", 1, KOKKOS_LAMBDA(const int) {
        corr_array_below(hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_eta) =
            corr_array_cloud(hm_metadata_iiPDF_Ncn, hm_metadata_iiPDF_eta);
      });
  // Estimates the correlation of the natural logarithm of a
  // hydrometeor species and eta using the correlation of chi and eta and the
  // correlation of chi and the natural logarithm of the hydrometeor.  This
  // facilitates the Cholesky decomposability of the correlation array that will
  // inevitably be decomposed for SILHS purposes. Without this estimation, we
  // have found that the resulting correlation matrix cannot be decomposed.
  for (jvar = ((hm_metadata_iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_9", 1, KOKKOS_LAMBDA(const int) {
          corr_array_cloud(jvar, hm_metadata_iiPDF_eta) =
              (corr_array_cloud(hm_metadata_iiPDF_eta, hm_metadata_iiPDF_chi) *
               corr_array_cloud(jvar, hm_metadata_iiPDF_chi));
        });
    Kokkos::parallel_for(
        "calc_corr_norm_and_cholesky_factor::scalar_10", 1, KOKKOS_LAMBDA(const int) {
          corr_array_below(jvar, hm_metadata_iiPDF_eta) =
              (corr_array_below(hm_metadata_iiPDF_eta, hm_metadata_iiPDF_chi) *
               corr_array_below(jvar, hm_metadata_iiPDF_chi));
        });
  }
  // Calc in cloud Cholesky
  Cholesky_factor(workspace, pdf_dim, corr_array_cloud,        // In
                  corr_array_scaling, corr_cholesky_mtx_cloud, // Out
                  l_corr_array_scaling                         // Out
  );
  // Calc out of cloud Cholesky
  Cholesky_factor(workspace, pdf_dim, corr_array_below,        // In
                  corr_array_scaling, corr_cholesky_mtx_below, // Out
                  l_corr_array_scaling                         // Out
  );
  // Use rc_1 to determine which correlation and Cholesky matrices to assign to 1st PDF
  // Correlation matrices are symmetric, so we copy ij values to ji indices as well
  // Cholesky matrices are lower triangular, so we only copy those values
  for (ivar = 0; ivar < pdf_dim; ++ivar) {
    for (jvar = ivar; jvar < pdf_dim; ++jvar) {
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "calc_corr_norm_and_cholesky_factor::cells_11", ((extent_0 * extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              if ((rc_1(i, k) > rc_tol)) {
                // Assign in cloud matrices to 1st PDF component
                corr_array_1_n(i, k, jvar, ivar)      = corr_array_cloud(jvar, ivar);
                corr_array_1_n(i, k, ivar, jvar)      = corr_array_cloud(jvar, ivar);
                corr_cholesky_mtx_1(i, k, jvar, ivar) = corr_cholesky_mtx_cloud(jvar, ivar);
              } else {
                // Assign out of cloud matrices to 1st PDF component
                corr_array_1_n(i, k, jvar, ivar)      = corr_array_below(jvar, ivar);
                corr_array_1_n(i, k, ivar, jvar)      = corr_array_below(jvar, ivar);
                corr_cholesky_mtx_1(i, k, jvar, ivar) = corr_cholesky_mtx_below(jvar, ivar);
              }
            });
      }
    }
  }
  // Use rc_2 to determine which correlation and Cholesky matrices to assign to 2nd PDF
  // Correlation matrices are symmetric, so we copy ij values to ji indices as well
  // Cholesky matrices are lower triangular, so we only copy those values
  for (ivar = 0; ivar < pdf_dim; ++ivar) {
    for (jvar = ivar; jvar < pdf_dim; ++jvar) {
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "calc_corr_norm_and_cholesky_factor::cells_12", ((extent_0 * extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              if ((rc_2(i, k) > rc_tol)) {
                // Assign in cloud matrices to 2nd PDF component
                corr_array_2_n(i, k, jvar, ivar)      = corr_array_cloud(jvar, ivar);
                corr_array_2_n(i, k, ivar, jvar)      = corr_array_cloud(jvar, ivar);
                corr_cholesky_mtx_2(i, k, jvar, ivar) = corr_cholesky_mtx_cloud(jvar, ivar);
              } else {
                // Assign out of cloud matrices to 2nd PDF component
                corr_array_2_n(i, k, jvar, ivar)      = corr_array_below(jvar, ivar);
                corr_array_2_n(i, k, ivar, jvar)      = corr_array_below(jvar, ivar);
                corr_cholesky_mtx_2(i, k, jvar, ivar) = corr_cholesky_mtx_below(jvar, ivar);
              }
            });
      }
    }
  }
  // Set upper triangular parts of Cholesky matrices to zero
  for (ivar = 0; ivar < pdf_dim; ++ivar) {
    for (jvar = 0; jvar < ivar; ++jvar) {
      {
        const int slice_extent_0 = corr_cholesky_mtx_1.extent_int(0);
        const int slice_extent_1 = corr_cholesky_mtx_1.extent_int(1);
        Kokkos::parallel_for(
            "calc_corr_norm_and_cholesky_factor::array_13", ((slice_extent_0 * slice_extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int slice_0                                 = ((cell % slice_extent_0));
              const int slice_1                                 = ((cell / slice_extent_0));
              corr_cholesky_mtx_1(slice_0, slice_1, jvar, ivar) = zero;
            });
      }
      {
        const int slice_extent_0 = corr_cholesky_mtx_2.extent_int(0);
        const int slice_extent_1 = corr_cholesky_mtx_2.extent_int(1);
        Kokkos::parallel_for(
            "calc_corr_norm_and_cholesky_factor::array_14", ((slice_extent_0 * slice_extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int slice_0                                 = ((cell % slice_extent_0));
              const int slice_1                                 = ((cell / slice_extent_0));
              corr_cholesky_mtx_2(slice_0, slice_1, jvar, ivar) = zero;
            });
      }
    }
  }
}
} // namespace
// =============================================================================
void comp_corr_norm(ScratchWorkspace& workspace,                                      // InOut
                    int nzt, int pdf_dim, int ngrdcol, int hydromet_dim,              // In
                    const hm_metadata_type& hm_metadata,                              // In
                    Real2D wm_zt, Real2D rc_1, Real2D rc_2,                           // In
                    Real2D mixt_frac,                                                 // In
                    Real2D precip_frac_1,                                             // In
                    Real2D precip_frac_2,                                             // In
                    Real2D wpNcnp_zt, Real3D wphydrometp_zt,                          // In
                    Real3D mu_x_1, Real3D mu_x_2, Real3D sigma_x_1, Real3D sigma_x_2, // In
                    Real3D sigma_x_1_n, Real3D sigma_x_2_n,                           // In
                    Real2D corr_array_n_cloud, Real2D corr_array_n_below,             // In
                    const pdf_parameter& pdf_params,                                  // In
                    int iiPDF_type,                                                   // In
                    bool l_calc_w_corr,                                               // In
                    bool l_fix_w_chi_eta_correlations,                                // In
                    Real4D corr_array_1_n, Real4D corr_array_2_n                      // Out
) {
  // Description:
  // Assemble symmetric normal-space correlation matrices from component
  // moments, prescribed cloud/clear values and the selected w-correlation
  // diagnosis. Retain each Native field-specific branch and source slot.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical levels
  // pdf_dim: Number of variables in the corr/mean/stdev arrays
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // wm_zt: Mean vertical velocity, <w>, on thermo. levels    [m/s]
  // rc_1: Mean of r_c (1st PDF component)                 [kg/kg]
  // rc_2: Mean of r_c (2nd PDF component)                 [kg/kg]
  // mixt_frac: Mixture fraction                                    [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)          [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)          [-]
  // wpNcnp_zt: Covariance of w and N_cn on t-levs.      [(m/s) num/kg]
  // wphydrometp_zt: Covariance of w and hm interp. to t-levs.  [(m/s)u.v.]
  // mu_x_1: Mean of x array (1st PDF component)          [units vary]
  // mu_x_2: Mean of x array (2nd PDF component)          [units vary]
  // sigma_x_1: Standard deviation of x array (1st PDF comp.)  [un. vary]
  // sigma_x_2: Standard deviation of x array (2nd PDF comp.)  [un. vary]
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // corr_array_n_cloud: Prescribed correlation array in cloud        [-]
  // corr_array_n_below: Prescribed correlation array below cloud     [-]
  // pdf_params: PDF parameters                               [units vary]
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF.
  // l_calc_w_corr: Calculate the correlations between w and the hydrometeors
  // l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta)
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]

  auto workspace_scope = workspace.scope();
  Real3D corr_w_hm_1_n = workspace_scope.view<Real3D>(
      "comp_corr_norm::corr_w_hm_1_n", ngrdcol, nzt,
      pdf_dim); // Correlation of w and ln hm (1st PDF component) ip   [-]
  Real3D corr_w_hm_2_n = workspace_scope.view<Real3D>(
      "comp_corr_norm::corr_w_hm_2_n", ngrdcol, nzt,
      pdf_dim); // Correlation of w and ln hm (2nd PDF component) ip   [-]
  Real2D corr_w_Ncn_1_n =
      workspace_scope.view<Real2D>("comp_corr_norm::corr_w_Ncn_1_n", ngrdcol,
                                   nzt); // Correlation of w and ln Ncn (1st PDF component)    [-]
  Real2D corr_w_Ncn_2_n =
      workspace_scope.view<Real2D>("comp_corr_norm::corr_w_Ncn_2_n", ngrdcol,
                                   nzt); // Correlation of w and ln Ncn (2nd PDF component)    [-]
  Real2D ones_vector =
      workspace_scope.view<Real2D>("comp_corr_norm::ones_vector", ngrdcol, nzt); // Vector of 1s
  Real2D Ncn_tol_in      = workspace_scope.view<Real2D>("comp_corr_norm::Ncn_tol_in", ngrdcol,
                                                        nzt); // Tolerance value for Ncn
  Real2D hydromet_tol_in = workspace_scope.view<Real2D>("comp_corr_norm::hydromet_tol_in", ngrdcol,
                                                        nzt); // Tolerance value for hydromet
  bool l_limit_corr_chi_eta{}; // Flag to limit the correlation of chi and eta [-]
  int ivar{};
  int jvar{};
  int hm_idx{};
  int i{};
  int k{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_w{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi     = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta     = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w       = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn     = hm_metadata.iiPDF_Ncn;
  const auto hm_metadata_hydromet_tol  = hm_metadata.hydromet_tol;
  const auto pdf_params_corr_chi_eta_1 = pdf_params.corr_chi_eta_1;
  const auto pdf_params_corr_chi_eta_2 = pdf_params.corr_chi_eta_2;
  const auto pdf_params_corr_w_chi_1   = pdf_params.corr_w_chi_1;
  const auto pdf_params_corr_w_chi_2   = pdf_params.corr_w_chi_2;
  const auto corr_array_n_below_host   = to_host_space(workspace, corr_array_n_below);
  const auto corr_array_n_cloud_host   = to_host_space(workspace, corr_array_n_cloud);

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // !! Normal space correlations
  // Initialize corr_w_hm_1_n and corr_w_hm_2_n arrays to 0.
  {
    const int slice_extent_0 = corr_w_hm_1_n.extent_int(0);
    const int slice_extent_1 = corr_w_hm_1_n.extent_int(1);
    const int slice_extent_2 = corr_w_hm_1_n.extent_int(2);
    Kokkos::parallel_for(
        "comp_corr_norm::array_1", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                        = ((cell % slice_extent_0));
          const int slice_1                        = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2                        = ((cell / (slice_extent_0 * slice_extent_1)));
          corr_w_hm_1_n(slice_0, slice_1, slice_2) = zero;
        });
  }
  {
    const int slice_extent_0 = corr_w_hm_2_n.extent_int(0);
    const int slice_extent_1 = corr_w_hm_2_n.extent_int(1);
    const int slice_extent_2 = corr_w_hm_2_n.extent_int(2);
    Kokkos::parallel_for(
        "comp_corr_norm::array_2", ((slice_extent_0 * slice_extent_1 * slice_extent_2)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                        = ((cell % slice_extent_0));
          const int slice_1                        = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2                        = ((cell / (slice_extent_0 * slice_extent_1)));
          corr_w_hm_2_n(slice_0, slice_1, slice_2) = zero;
        });
  }
  // Set ones_vector to a vector of 1s.
  {
    const int slice_extent_0 = ones_vector.extent_int(0);
    const int slice_extent_1 = ones_vector.extent_int(1);
    Kokkos::parallel_for(
        "comp_corr_norm::array_3", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0             = ((cell % slice_extent_0));
          const int slice_1             = ((cell / slice_extent_0));
          ones_vector(slice_0, slice_1) = one;
        });
  }
  // Calculate normal space correlations involving w by first calculating total
  // covariances involving w (<w'Ncn'>, etc.) using the down-gradient
  // approximation.
  if (l_calc_w_corr) {
    {
      const int slice_extent_0 = Ncn_tol_in.extent_int(0);
      const int slice_extent_1 = Ncn_tol_in.extent_int(1);
      Kokkos::parallel_for(
          "comp_corr_norm::array_4", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0            = ((cell % slice_extent_0));
            const int slice_1            = ((cell / slice_extent_0));
            Ncn_tol_in(slice_0, slice_1) = Ncn_tol;
          });
    }
    // Calculate the correlation of w and ln Ncn in each PDF component.
    // The subroutine calc_corr_w_hm_n can be used to do this as long as a
    // value of 1 is sent in for precip_frac_1 and precip_frac_2.
    jvar = iiPDF_Ncn;
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "comp_corr_norm::cells_5", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            calc_corr_w_hm_n(wm_zt(i, k), wpNcnp_zt(i, k),              // In
                             mu_x_1(i, k, iiPDF_w),                     // In
                             mu_x_2(i, k, iiPDF_w),                     // In
                             mu_x_1(i, k, jvar), mu_x_2(i, k, jvar),    // In
                             sigma_x_1(i, k, iiPDF_w),                  // In
                             sigma_x_2(i, k, iiPDF_w),                  // In
                             sigma_x_1(i, k, jvar),                     // In
                             sigma_x_2(i, k, jvar),                     // In
                             sigma_x_1_n(i, k, jvar),                   // In
                             sigma_x_2_n(i, k, jvar),                   // In
                             mixt_frac(i, k), ones_vector(i, k),        // In
                             ones_vector(i, k), Ncn_tol_in(i, k),       // In
                             corr_w_Ncn_1_n(i, k), corr_w_Ncn_2_n(i, k) // Out
            );
          });
    }
    // Calculate the correlation of w and the natural logarithm of the
    // hydrometeor for each PDF component and each hydrometeor type.
    for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
      hm_idx = (pdf2hydromet_idx(jvar, hm_metadata) + 1) - 1;
      {
        const int slice_extent_0 = hydromet_tol_in.extent_int(0);
        const int slice_extent_1 = hydromet_tol_in.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_6", ((slice_extent_0 * slice_extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int slice_0                 = ((cell % slice_extent_0));
              const int slice_1                 = ((cell / slice_extent_0));
              hydromet_tol_in(slice_0, slice_1) = hm_metadata_hydromet_tol(hm_idx);
            });
      }
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "comp_corr_norm::cells_7", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              calc_corr_w_hm_n(wm_zt(i, k), wphydrometp_zt(i, k, hm_idx),  // In
                               mu_x_1(i, k, iiPDF_w),                      // In
                               mu_x_2(i, k, iiPDF_w),                      // In
                               mu_x_1(i, k, jvar), mu_x_2(i, k, jvar),     // In
                               sigma_x_1(i, k, iiPDF_w),                   // In
                               sigma_x_2(i, k, iiPDF_w),                   // In
                               sigma_x_1(i, k, jvar),                      // In
                               sigma_x_2(i, k, jvar),                      // In
                               sigma_x_1_n(i, k, jvar),                    // In
                               sigma_x_2_n(i, k, jvar),                    // In
                               mixt_frac(i, k), precip_frac_1(i, k),       // In
                               precip_frac_2(i, k), hydromet_tol_in(i, k), // In
                               corr_w_hm_1_n(i, k, jvar),                  // Out
                               corr_w_hm_2_n(i, k, jvar)                   // Out
              );
            });
      }
    }
    // Native loop: jvar = iiPDF_Ncn+1, pdf_dim
  }
  // l_calc_w_corr
  // In order to decompose the normal space correlation matrix,
  // we must not have a perfect correlation of chi and
  // eta. Thus, we impose a limitation.
  l_limit_corr_chi_eta = true;
  // Initialize the normal space correlation arrays
  {
    const int slice_extent_0 = corr_array_1_n.extent_int(0);
    const int slice_extent_1 = corr_array_1_n.extent_int(1);
    const int slice_extent_2 = corr_array_1_n.extent_int(2);
    const int slice_extent_3 = corr_array_1_n.extent_int(3);
    Kokkos::parallel_for(
        "comp_corr_norm::array_8",
        ((slice_extent_0 * slice_extent_1 * slice_extent_2 * slice_extent_3)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2 = (((cell / (slice_extent_0 * slice_extent_1)) % slice_extent_2));
          const int slice_3 = ((cell / (slice_extent_0 * slice_extent_1 * slice_extent_2)));
          corr_array_1_n(slice_0, slice_1, slice_2, slice_3) = zero;
        });
  }
  {
    const int slice_extent_0 = corr_array_2_n.extent_int(0);
    const int slice_extent_1 = corr_array_2_n.extent_int(1);
    const int slice_extent_2 = corr_array_2_n.extent_int(2);
    const int slice_extent_3 = corr_array_2_n.extent_int(3);
    Kokkos::parallel_for(
        "comp_corr_norm::array_9",
        ((slice_extent_0 * slice_extent_1 * slice_extent_2 * slice_extent_3)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = (((cell / slice_extent_0) % slice_extent_1));
          const int slice_2 = (((cell / (slice_extent_0 * slice_extent_1)) % slice_extent_2));
          const int slice_3 = ((cell / (slice_extent_0 * slice_extent_1 * slice_extent_2)));
          corr_array_2_n(slice_0, slice_1, slice_2, slice_3) = zero;
        });
  }
  // !! The corr_arrays are assumed to be lower triangular matrices
  // Set diagonal elements to 1
  for (ivar = 0; ivar < pdf_dim; ++ivar) {
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "comp_corr_norm::cells_10", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
            const int i                      = ((cell % extent_1));
            const int k                      = ((cell / extent_1));
            corr_array_1_n(i, k, ivar, ivar) = one;
            corr_array_2_n(i, k, ivar, ivar) = one;
          });
    }
  }
  // !! This code assumes the following order in the prescribed correlation
  // !! arrays (iiPDF indices):
  // !! chi, eta, w, Ncn, <hydrometeors> (indices increasing from left to right)
  if (l_fix_w_chi_eta_correlations) {
    // Correlation of chi (old s) and eta (old t)
    component_corr_chi_eta(nzt, ngrdcol,                                  // In
                           rc_1,                                          // In
                           rc_2,                                          // In
                           corr_array_n_cloud_host(iiPDF_eta, iiPDF_chi), // In
                           corr_array_n_below_host(iiPDF_eta, iiPDF_chi), // In
                           l_limit_corr_chi_eta,                          // In
                           Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_eta,
                                           iiPDF_chi), // Out
                           Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_eta,
                                           iiPDF_chi) // Out
    );
  } else {
    // Preferred, more accurate version.
    for (i = 0; i < ngrdcol; ++i) {
      {
        const int slice_extent_0 = corr_array_1_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_11", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0 = cell;
              corr_array_1_n(i, slice_0, iiPDF_eta, iiPDF_chi) =
                  pdf_params_corr_chi_eta_1(i, slice_0);
            });
      }
      {
        const int slice_extent_0 = corr_array_2_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_12", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0 = cell;
              corr_array_2_n(i, slice_0, iiPDF_eta, iiPDF_chi) =
                  pdf_params_corr_chi_eta_2(i, slice_0);
            });
      }
    }
  }
  if (l_fix_w_chi_eta_correlations) {
    // Correlation of chi (old s) and w
    component_corr_w_x(
        nzt, ngrdcol,                                // In
        rc_1,                                        // In
        rc_2,                                        // In
        corr_array_n_cloud_host(iiPDF_w, iiPDF_chi), // In
        corr_array_n_below_host(iiPDF_w, iiPDF_chi), // In
        iiPDF_type,                                  // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_w,
                        iiPDF_chi),                                                   // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_chi) // Out
    );
  } else {
    // Preferred, more accurate version.
    for (i = 0; i < ngrdcol; ++i) {
      {
        const int slice_extent_0 = corr_array_1_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_13", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0                              = cell;
              corr_array_1_n(i, slice_0, iiPDF_w, iiPDF_chi) = pdf_params_corr_w_chi_1(i, slice_0);
            });
      }
      {
        const int slice_extent_0 = corr_array_2_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_14", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0                              = cell;
              corr_array_2_n(i, slice_0, iiPDF_w, iiPDF_chi) = pdf_params_corr_w_chi_2(i, slice_0);
            });
      }
    }
  }
  // Correlation of chi (old s) and ln Ncn, corr_array_n_cloud used twice because
  // Ncn is an inherently in-cloud property.
  component_corr_x_hm_n_ip(
      nzt, ngrdcol,                                  // In
      rc_1,                                          // In
      rc_2,                                          // In
      corr_array_n_cloud_host(iiPDF_Ncn, iiPDF_chi), // In
      corr_array_n_cloud_host(iiPDF_Ncn, iiPDF_chi), // In
      Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                      iiPDF_chi),                                                     // Out
      Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi) // Out
  );
  // Correlation of chi (old s) and the natural logarithm of the hydrometeors
  for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    component_corr_x_hm_n_ip(
        nzt, ngrdcol,                                                               // In
        rc_1,                                                                       // In
        rc_2,                                                                       // In
        corr_array_n_cloud_host(jvar, iiPDF_chi),                                   // In
        corr_array_n_below_host(jvar, iiPDF_chi),                                   // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_chi), // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_chi)  // Out
    );
  }
  // Correlation of eta (old t) and w
  if (l_fix_w_chi_eta_correlations) {
    // Correlation of chi (old s) and w
    component_corr_w_x(
        nzt, ngrdcol,                                // In
        rc_1,                                        // In
        rc_2,                                        // In
        corr_array_n_cloud_host(iiPDF_w, iiPDF_eta), // In
        corr_array_n_below_host(iiPDF_w, iiPDF_eta), // In
        iiPDF_type,                                  // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_w,
                        iiPDF_eta),                                                   // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_eta) // Out
    );
  } else {
    // Preferred, more accurate version.
    for (i = 0; i < ngrdcol; ++i) {
      {
        const int slice_extent_0 = corr_array_1_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_15", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0                              = cell;
              corr_array_1_n(i, slice_0, iiPDF_w, iiPDF_chi) = pdf_params_corr_w_chi_1(i, slice_0);
            });
      }
      {
        const int slice_extent_0 = corr_array_2_n.extent_int(1);
        Kokkos::parallel_for(
            "comp_corr_norm::array_16", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
              const int slice_0                              = cell;
              corr_array_2_n(i, slice_0, iiPDF_w, iiPDF_chi) = pdf_params_corr_w_chi_2(i, slice_0);
            });
      }
    }
  }
  // Correlation of eta (old t) and ln Ncn, corr_array_n_cloud used twice because
  // Ncn is an inherently in-cloud property.
  component_corr_x_hm_n_ip(
      nzt, ngrdcol,                                  // In
      rc_1,                                          // In
      rc_2,                                          // In
      corr_array_n_cloud_host(iiPDF_Ncn, iiPDF_eta), // In
      corr_array_n_cloud_host(iiPDF_Ncn, iiPDF_eta), // In
      Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                      iiPDF_eta),                                                     // Out
      Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta) // Out
  );
  // Correlation of eta (old t) and the natural logarithm of the hydrometeors
  for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    component_corr_eta_hm_n_ip(
        nzt, ngrdcol, // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_eta,
                        iiPDF_chi),                                                 // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_chi), // In
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_eta,
                        iiPDF_chi),                                                 // In
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_chi), // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_eta), // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_eta)  // Out
    );
  }
  // Correlation of w and ln Ncn
  component_corr_w_hm_n_ip(
      nzt, ngrdcol,                                                                  // In
      corr_w_Ncn_1_n, rc_1,                                                          // In
      corr_w_Ncn_2_n, rc_2,                                                          // In
      corr_array_n_cloud_host(iiPDF_Ncn, iiPDF_w),                                   // In
      corr_array_n_below_host(iiPDF_Ncn, iiPDF_w),                                   // In
      l_calc_w_corr,                                                                 // In
      Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // Out
      Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w)  // Out
  );
  // Correlation of w and the natural logarithm of the hydrometeors
  for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    component_corr_w_hm_n_ip(
        nzt, ngrdcol,                                                             // In
        Kokkos::subview(corr_w_hm_1_n, Kokkos::ALL, Kokkos::ALL, jvar), rc_1,     // In
        Kokkos::subview(corr_w_hm_2_n, Kokkos::ALL, Kokkos::ALL, jvar), rc_2,     // In
        corr_array_n_cloud_host(jvar, iiPDF_w),                                   // In
        corr_array_n_below_host(jvar, iiPDF_w),                                   // In
        l_calc_w_corr,                                                            // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_w), // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_w)  // Out
    );
  }
  // Correlation of ln Ncn and the natural logarithm of the hydrometeors
  for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    component_corr_hmx_hmy_n_ip(
        nzt, ngrdcol,                                                               // In
        rc_1,                                                                       // In
        rc_2,                                                                       // In
        corr_array_n_cloud_host(jvar, iiPDF_Ncn),                                   // In
        corr_array_n_below_host(jvar, iiPDF_Ncn),                                   // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn), // Out
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn)  // Out
    );
  }
  // Correlation of the natural logarithm of two hydrometeors
  for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim - 1; ++ivar) {
    for (jvar = ((ivar + 1)); jvar < pdf_dim; ++jvar) {
      component_corr_hmx_hmy_n_ip(
          nzt, ngrdcol,                                                          // In
          rc_1,                                                                  // In
          rc_2,                                                                  // In
          corr_array_n_cloud_host(jvar, ivar),                                   // In
          corr_array_n_below_host(jvar, ivar),                                   // In
          Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // Out
          Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar)  // Out
      );
    }
    // jvar
  }
  // ivar
  // For ease of use later in the code, we make the correlation arrays
  // symmetrical
  for (k = 0; k < nzt; ++k) {
    for (i = 0; i < ngrdcol; ++i) {
      mirror_lower_triangular_matrix(
          pdf_dim,                                                        // In
          Kokkos::subview(corr_array_1_n, i, k, Kokkos::ALL, Kokkos::ALL) // InOut
      );
      mirror_lower_triangular_matrix(
          pdf_dim,                                                        // In
          Kokkos::subview(corr_array_2_n, i, k, Kokkos::ALL, Kokkos::ALL) // InOut
      );
    }
  }
  return;
}
// =============================================================================
void calc_comp_mu_sigma_hm(int nzt, int ngrdcol,                       // In
                           Real2D hmm, Real2D hmp2,                    // In
                           real_t hmp2_ip_on_hmm2_ip,                  // In
                           Real2D mixt_frac, Real2D precip_frac,       // In
                           Real2D precip_frac_1, Real2D precip_frac_2, // In
                           real_t hm_tol, Real1D precip_frac_tol,      // In
                           Real2D mu_thl_1, Real2D mu_thl_2,           // In
                           Real1D omicron, Real1D zeta_vrnce_rat_in,   // In
                           Real2D mu_hm_1, Real2D mu_hm_2,             // Out
                           Real2D sigma_hm_1, Real2D sigma_hm_2,       // Out
                           Real2D hm_1, Real2D hm_2,                   // Out
                           Real2D sigma_hm_1_sqd_on_mu_hm_1_sqd,       // Out
                           Real2D sigma_hm_2_sqd_on_mu_hm_2_sqd        // Out
) {
  // Description:
  // Calculates the in-precipitation mean and in-precipitation standard
  // deviation in both PDF components for a precipitating hydrometeor.
  // When precipitation is found in both PDF components (precip_frac_1 > 0 and
  // precip_frac_2 > 0), the method that solves for in-precip. mean and
  // in-precip. standard deviation in each PDF component, preserving overall
  // mean and overall variance, is used.  When precipitation fraction is found
  // in one PDF component but not the other one (precip_frac_1 > 0 and
  // precip_frac_2 = 0, or precip_frac_1 = 0 and precip_frac_2 > 0), the
  // calculation of component in-precip. mean and in-precip. standard deviation
  // is simple.  When precipitation is not found in either component
  // (precip_frac_1 = 0 and precip_frac_2 = 0), there isn't any precipitation
  // found overall (at that grid level).
  // DESCRIPTION OF THE METHOD THAT SOLVES FOR TWO IN-PRECIPITATION COMPONENTS
  // =========================================================================
  // OVERVIEW
  // The goal is to calculate the in-precip. mean of the hydrometeor field in
  // each PDF component (mu_hm_1 and mu_hm_2) in a scenario when there is
  // precipitation found in both PDF components.  The fields provided are the
  // overall mean of the hydrometeor, <hm>, the overall variance of the
  // hydrometeor, <hm’^2>, the mixture fraction, a, the overall precipitation
  // fraction, f_p, and the precipitation fraction in each PDF component
  // (f_p_1 and f_p_2).
  // The PDF equation for <hm> is:
  // <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2.
  // Likewise, the PDF equation for <hm’^2> is:
  // <hm’^2> = a * f_p_1 * ( mu_hm_1^2 + sigma_hm_1^2 )
  //           + ( 1 - a ) * f_p_2 * ( mu_hm_2^2 + sigma_hm_2^2 )
  //           - <hm>^2;
  // where sigma_hm_1 and sigma_hm_2 are the in-precip. standard deviations of
  // the hydrometeor field in each PDF component.  This can be rewritten as:
  // <hm’^2>
  // = a * f_p_1 * ( 1 + sigma_hm_1^2 / mu_hm_1^2 ) * mu_hm_1^2
  //   + ( 1 - a ) * f_p_2 * ( 1 + sigma_hm_2^2 / mu_hm_2^2 ) * mu_hm_2^2
  //   - <hm>^2.
  // The ratio of sigma_hm_2^2 to mu_hm_2^2 is denoted R:
  // R = sigma_hm_2^2 / mu_hm_2^2.
  // In order to allow sigma_hm_1^2 / mu_hm_1^2 to have a different ratio, the
  // parameter zeta is introduced, such that:
  // R * ( 1 + zeta ) = sigma_hm_1^2 / mu_hm_1^2;
  // where zeta > -1.  When -1 < zeta < 0, the ratio sigma_hm_2^2 / mu_hm_2^2
  // grows at the expense of sigma_hm_1^2 / mu_hm_1^2, which narrows.  When
  // zeta = 0, the ratio sigma_hm_1^2 / mu_hm_1^2 is the same as
  // sigma_hm_2^2 / mu_hm_2^2.  When zeta > 0, sigma_hm_1^2 / mu_hm_1^2 grows
  // at the expense of sigma_hm_2^2 / mu_hm_2^2, which narrows.  The component
  // variances are written as:
  // sigma_hm_1^2 = R * ( 1 + zeta ) * mu_hm_1^2; and
  // sigma_hm_2^2 = R * mu_hm_2^2,
  // and the component standard deviations are simply:
  // sigma_hm_1 = sqrt( R * ( 1 + zeta ) ) * mu_hm_1; and
  // sigma_hm_2 = sqrt( R ) * mu_hm_2.
  // The equation for <hm’^2> can be rewritten as:
  // <hm’^2> = a * f_p_1 * ( 1 + R * ( 1 + zeta ) ) * mu_hm_1^2
  //           + ( 1 - a ) * f_p_2 * ( 1 + R ) * mu_hm_2^2
  //           - <hm>^2.
  // HYDROMETEOR IN-PRECIP. VARIANCE:
  // THE SPREAD OF THE MEANS VS. THE STANDARD DEVIATIONS
  // Part I:  Minimum and Maximum Values for R
  // The in-precip. variance of the hydrometeor is accounted for through a
  // combination of the variance of each PDF component and the spread between
  // the means of each PDF component.  At one extreme, the standard deviation
  // of each component could be set to 0 and the in-prccip. variance could be
  // accounted for by spreading the PDF component (in-precip.) means far apart.
  // The value of R in this scenario would be its minimum possible value, which
  // is 0.  At the other extreme, the means of each component could be set
  // equal to each other and the in-precip. variance could be accounted for
  // entirely by the PDF component (in-precip.) standard deviations.  The value
  // of R in this scenario would be its maximum possible value, which is Rmax.
  // In order to calculate the value of Rmax, use the equation set but set
  // mu_hm_1 = mu_hm_2 and R = Rmax.  When this happens:
  // <hm> = ( a * f_p_1 + ( 1- a ) * f_p_2 ) * mu_hm_i;
  // and since f_p = a * f_p_1 + ( 1 - a ) * f_p_2:
  // mu_hm_i = <hm> / f_p = <hm|_ip>;
  // where <hm|_ip> is the in-precip. mean of the hydrometeor.  The equation
  // for hydrometeor variance in this scenario becomes:
  // <hm’^2> = <hm|_ip>^2 * ( a * f_p_1 * ( 1 + Rmax * ( 1 + zeta ) )
  //                          + ( 1 - a ) * f_p_2 * ( 1 + Rmax ) )
  //           - <hm>^2.
  // The general equation for the in-precip. variance of a hydrometeor,
  // <hm|_ip’^2>, is given by:
  // <hm|_ip’^2> = ( <hm’^2> + <hm>^2 - f_p * <hm|_ip>^2 ) / f_p;
  // which can be rewritten as:
  // <hm’^2> + <hm>^2 = f_p * ( <hm|_ip’^2> + <hm|_ip>^2 ).
  // When the above equation is substituted into the modified PDF equation for
  // <hm’^2>, Rmax is solved for and the equation is:
  // Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
  //        * ( <hm|_ip’^2> / <hm|_ip>^2 ).
  // Here, in the scenario that zeta = 0, both PDF components have the same
  // mean and same variance, which reduces the in-precip. distribution to an
  // assumed single lognormal, and the above equation reduces to:
  // Rmax = <hm|_ip’^2> / <hm|_ip>^2;
  // which is what is expected in that case.
  // Part II:  Enter omicron
  // A parameter is used to prescribe the ratio of R to its maximum value,
  // Rmax.  The prescribed parameter is called omicron, where:
  // R = omicron * Rmax;
  // where 0 <= omicron <= 1.  When omicron = 0, the standard deviation of each
  // PDF component is 0, and mu_hm_1 is spread as far away from mu_hm_2 as it
  // needs to be to account for the in-precip. variance.  When omicron = 1,
  // mu_hm_1 is equal to mu_hm_2, and the standard deviations of the PDF
  // components account for all of the in-precip. variance (and when zeta = 0,
  // the PDF shape is a single lognormal in-precip.).  At intermediate values
  // of omicron, the means of each PDF component are somewhat spread and each
  // PDF component has some width.  The modified parameters are listed below.
  // The ratio of sigma_hm_2^2 to mu_hm_2^2 is:
  // sigma_hm_2^2 / mu_hm_2^2 = omicron * Rmax;
  // and the ratio of sigma_hm_1^2 / mu_hm_1^2 is:
  // sigma_hm_1^2 / mu_hm_1^2 = omicron * Rmax * ( 1 + zeta ).
  //
  // The component variances are written as:
  // sigma_hm_1^2 = omicron * Rmax * ( 1 + zeta ) * mu_hm_1^2; and
  // sigma_hm_2^2 = omicron * Rmax * mu_hm_2^2,
  // and the component standard deviations are simply:
  // sigma_hm_1 = sqrt( omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1; and
  // sigma_hm_2 = sqrt( omicron * Rmax ) * mu_hm_2.
  // The equation set becomes:
  // [1] <hm> = a * f_p_1 * mu_hm_1 + ( 1- a ) * f_p_2 * mu_hm_2; and
  // [2] <hm’^2>
  //     = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) ) * mu_hm_1^2
  //       + ( 1 - a ) * f_p_2 * ( 1 + omicron * Rmax ) * mu_hm_2^2
  //       - <hm>^2.
  // SOLVING THE EQUATION SET FOR MU_HM_1 AND MU_HM_2.
  // The above system of two equations can be solved for mu_hm_1 and mu_hm_2.
  // All other quantities in the equation set are known quantities.  The
  // equation for <hm> is rewritten to isolate mu_hm_2:
  // mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).
  // The above equation is substituted into the equation for <hm’^2>.  The
  // equation for <hm’^2> is rewritten, resulting in:
  // [ a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
  //   + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ]
  // * mu_hm_1^2
  // + [ - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
  //     / ( ( 1 - a ) * f_p_2 ) ] * mu_hm_1
  // + [ - ( <hm’^2>
  //         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
  //           * <hm>^2 ) ]
  // = 0.
  // This equation is of the form:
  // A * mu_hm_1^2 + B * mu_hm_1 + C = 0;
  // so the solution for mu_hm_1 is:
  // mu_hm_1 = ( -B +/- sqrt( B^2 - 4*A*C ) ) / (2*A);
  // where:
  // A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
  //     + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 );
  // B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
  //     / ( ( 1 - a ) * f_p_2 );
  // and
  // C = - ( <hm’^2>
  //         + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
  //           * <hm>^2 ).
  // The signs of the coefficients:
  // 1) coefficient A is always positive,
  // 2) coefficient B is always negative (this means that -B is always
  //    positive), and
  // 3) coefficient C can be positive, negative, or zero.
  // Since ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) ) * <hm>^2 is
  // always negative and <hm’^2> is always positive, the sign of coefficient C
  // depends on which term is greater in magnitude.
  // When <hm’^2> is greater, the sign of coefficient C is negative.  This
  // means that -4*A*C is positive, which in turn means that
  // sqrt( B^2 - 4*A*C ) is greater in magnitude than -B.  If the subtraction
  // option of the +/- were to be chosen, the value of mu_hm_1 would be
  // negative in this scenerio.  So the natural thing to do would be to always
  // choose the addition option.  However, this method requires that mu_hm_1
  // equals mu_hm_2 when omicron = 1.  When zeta >= 0, this happens when the
  // addition option is chosen, but not when the subtraction option is chosen.
  // However, when zeta < 0, this happens when the subtraction option is
  // chosen, but not when the addition option is chosen.  So, the equation for
  // mu_hm_1 becomes:
  // mu_hm_1 = ( -B + sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta >= 0; and
  // mu_hm_1 = ( -B - sqrt( B^2 - 4*A*C ) ) / (2*A); when zeta < 0.
  // Once this is set, of course:
  // mu_hm_2 = ( <hm> - a * f_p_1 * mu_hm_1 ) / ( ( 1 - a ) * f_p_2 ).
  // The system has been solved and the in-precip. PDF component means have
  // been found!
  // NOTES
  // Note 1:
  // The term B^2 - 4*A*C has been analyzed, and mathematically:
  // B^2 - 4*A*C >= 0
  // always holds true.  Additionally, the minimum value:
  // B^2 - 4*A*C = 0,
  // can only occur when omicron = 1 and zeta = 0 (or alternatively to
  // zeta = 0, Rmax = 0, but this only occurs when <hm|_ip'^2> / <hm|_ip>^2 has
  // a value of 0).
  // Numerically, when omicron = 1 and zeta = 0, B^2 - 4*A*C can produce very
  // small (on the order of epsilon) negative values.  This is due to numerical
  // round off error.  When this happens, the erroneous small, negative value
  // of B^2 - 4*A*C is simply reset to the value it's supposed to have, which
  // is 0.
  // Note 2:
  // As the value of <hm|_ip'^2> / <hm|_ip>^2 increases and as the value of
  // omicron decreases (narrowing the in-precip standard deviations and
  // increasing the spread between the in-precip means), a situtation arises
  // where the value of one of the component means will become negative.  This
  // is because there is a limit to the amount of in-precip variance that can
  // be represented by this kind of distribution.  In order to prevent
  // out-of-bounds values of mu_hm_1 or mu_hm_2, lower limits will be
  // declared, called mu_hm_1_min and mu_hm_2_min.  The value of the
  // hydrometeor in-precip. component mean will be limited from going any
  // smaller (or negative) at this value.  From there, the value of the other
  // hydrometeor in-precip. component mean is easy to calculate.  Then, both
  // values will be entered into the calculation of hydrometeor variance, which
  // will be rewritten to solve for R.  Then, both the hydrometeor mean and
  // hydrometeor variance will be preserved with a valid distribution.
  // In this emergency scenario, the value of R is:
  // R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
  //       - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
  //     / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
  //         + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
  // The minimum values of the in-precip. component means are bounded by:
  // mu_hm_1_min >= hm_tol / f_p_1; and
  // mu_hm_2_min >= hm_tol / f_p_2.
  // These are set this way because hm_1 ( = mu_hm_1 * f_p_1 ) and
  // hm_2 ( = mu_hm_2 * f_p_2 ) need to have values of at least hm_tol when
  // precipitation is found in both PDF components.
  // However, an in-precip. component mean value of hm_tol / f_p_1 or
  // hm_tol / f_p_2 often produces a distribution where one component centers
  // around values that are too small to be a good match with data taken from
  // Large Eddy Simulations (LES).  It is desirable to increase the minimum
  // threshold of mu_hm_1 and mu_hm_2.
  // As the minimum threshold increases, the value of the in-precip. component
  // mean that is from the component that is not being set to the minimum
  // threshold decreases.  If the minimum threshold were to be boosted as high
  // as <hm> / f_p (in most cases, <hm> / f_p >> hm_tol / f_p_i), both
  // components would have a value of <hm> / f_p.  The minimum threshold should
  // not be set this high.
  // Additionally, the minimum threshold for one in-precip. component mean
  // cannot be set so high as to drive the other in-precip. component mean
  // below hm_tol / f_p_i.  (This doesn't come into play unless <hm> is close
  // to hm_tol.)  The upper limit for the in-precip. mean values are:
  // mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * f_p_2 * ( hm_tol / f_p_2 ) )
  //                          / ( a * f_p_1 ); and
  // mu_hm_2|_(upper. lim.) = ( <hm> - a * f_p_1 * ( hm_tol / f_p_1 ) )
  //                          / ( ( 1 - a ) * f_p_2 );
  // which reduces to:
  // mu_hm_1|_(upper. lim.) = ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 );
  // and
  // mu_hm_2|_(upper. lim.) = ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ).
  // An appropriate minimum value for mu_hm_1 can be set by:
  // mu_hm_1_min = | min( hm_tol / f_p_1
  //               |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_1 ),
  //               |      ( <hm> - ( 1 - a ) * hm_tol ) / ( a * f_p_1 ) );
  //               |    where <hm> / f_p > hm_tol / f_p_1;
  //               | hm_tol / f_p_1;
  //               |    where <hm> / f_p <= hm_tol / f_p_1;
  // and similarly for mu_hm_2:
  // mu_hm_2_min = | min( hm_tol / f_p_2
  //               |      + mu_hm_min_coef * ( <hm> / f_p - hm_tol / f_p_2 ),
  //               |      ( <hm> - a * hm_tol ) / ( ( 1 - a ) * f_p_2 ) );
  //               |    where <hm> / f_p > hm_tol / f_p_2;
  //               | hm_tol / f_p_2;
  //               |    where <hm> / f_p <= hm_tol / f_p_2;
  // where mu_hm_min_coef is a coefficient that has a value
  // 0 <= mu_hm_min_coef < 1.  When the value of mu_hm_min_coef is 0,
  // mu_hm_1_min reverts to hm_tol / f_p_1 and mu_hm_2_min reverts to
  // hm_tol / f_p_2.  An appropriate value for mu_hm_min_coef should be small,
  // such as 0.01 - 0.05.
  // Note 3:
  // When the value of zeta >= 0, the value of mu_hm_1 tends to be larger than
  // the value of mu_hm_2.  Likewise when the value of zeta < 0, the value of
  // mu_hm_2 tends to be larger than the value of mu_hm_1.  Since most cloud
  // water and cloud fraction tends to be found in PDF component 1, it is
  // advantageous to have the larger in-precip. component mean of the
  // hydrometeor also found in PDF component 1.  The recommended value of zeta
  // is a value greater than or equal to 0.
  // Update:
  // In order to better represent the increase in <th_l'^2> near the ground
  // from the evaporation of rain, the code is modified to tend towards a
  // negative correlation of th_l and hm.  When mu_thl_1 <= mu_thl_2,
  // mu_hm_1 >= mu_hm_2; otherwise, when mu_thl_1 > mu_thl_2,
  // mu_hm_1 <= mu_hm_2.
  // In the original derivation, mu_hm_1 >= mu_hm_2 when zeta >= 0, and
  // mu_hm_1 < mu_hm_2 when zeta < 0, where zeta is a tunable or adjustable
  // parameter.  In order to allow the relationship of mu_hm_1 to mu_hm_2 to
  // depend on the relationship of mu_thl_1 to mu_thl_2, the value of zeta must
  // also depend on the relationship of mu_thl_1 to mu_thl_2.
  // When mu_thl_1 <= mu_thl_2:
  // The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 >= mu_hm_2, so
  // zeta >= 0.  The tunable value of zeta is referred to as zeta_in.  When
  // zeta_in is already greater than 0 (meaning sigma_hm_1^2 / mu_hm_1^2 is
  // greater than sigma_hm_2^2 / mu_hm_2^2), zeta is simply set to zeta_in.  In
  // other words, the component with the greater mean value of the hydrometeor
  // also has the greater value of component variance to the square of the
  // component mean.  However, when zeta_in is less than 0, zeta must be
  // adjusted to be greater than 0.  The following equation is used to set the
  // value of zeta:
  //        | zeta_in; when zeta_in >= 0
  // zeta = |
  //        | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in < 0.
  // Previously, when zeta_in < 0, mu_hm_1 < mu_hm_2, and
  // sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
  // be greater than 0, sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.
  // The ratio of the greater variance-over-mean-squared to the smaller
  // variance-over-mean-squared remains the same when using the equation for
  // zeta listed above.
  // When mu_thl_1 > mu_thl_2:
  // The relationship of mu_hm_1 to mu_hm_2 is mu_hm_1 <= mu_hm_2, so
  // zeta <= 0.  When zeta_in is already less than 0 (meaning
  // sigma_hm_1^2 / mu_hm_1^2 is less than sigma_hm_2^2 / mu_hm_2^2), zeta is
  // simply set to zeta_in.  In other words, the component with the greater
  // mean value of the hydrometeor also has the greater value of component
  // variance to the square of the component mean.  However, when zeta_in is
  // greater than 0, zeta must be adjusted to be less than 0.  The following
  // equation is used to set the value of zeta:
  //        | zeta_in; when zeta_in <= 0
  // zeta = |
  //        | ( 1 / ( 1 + zeta_in ) ) - 1; when zeta_in > 0.
  // Previously, when zeta_in > 0, mu_hm_1 > mu_hm_2, and
  // sigma_hm_1^2 / mu_hm_1^2 > sigma_hm_2^2 / mu_hm_2^2.  Now that zeta has to
  // be less than 0, sigma_hm_1^2 / mu_hm_1^2 < sigma_hm_2^2 / mu_hm_2^2.
  // The ratio of the greater variance-over-mean-squared to the smaller
  // variance-over-mean-squared remains the same when using the equation for
  // zeta listed above.
  // Brian Griffin; February 2015.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hmm: Hydrometeor mean (overall), <hm>           [hm units]
  // hmp2: Hydrometeor variance (overall), <hm'^2>    [hm units^2]
  // hmp2_ip_on_hmm2_ip: Ratio <hm|_ip'^2> / <hm|_ip>^2                 [-]
  // mixt_frac: Mixture fraction                           [-]
  // precip_frac: Precipitation fraction (overall)           [-]
  // precip_frac_1: Precipitation fraction (1st PDF component) [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component) [-]
  // hm_tol: Tolerance value of hydrometeor          [hm units]
  // precip_frac_tol: Min. precip. frac. when hydromet. are present  [-]
  // mu_thl_1: Mean of th_l (1st PDF component)           [K]
  // mu_thl_2: Mean of th_l (2nd PDF component)           [K]
  // omicron: Relative width parameter, omicron = R / Rmax    [-]
  // zeta_vrnce_rat_in: Width parameter for sigma_hm_1^2 / mu_hm_1^2    [-]
  // mu_hm_1: Mean of hm (1st PDF component) in-precip (ip)      [hm un]
  // mu_hm_2: Mean of hm (2nd PDF component) ip                  [hm un]
  // sigma_hm_1: Standard deviation of hm (1st PDF component) ip    [hm un]
  // sigma_hm_2: Standard deviation of hm (2nd PDF component) ip    [hm un]
  // hm_1: Mean of hm (1st PDF component)                     [hm un]
  // hm_2: Mean of hm (2nd PDF component)                     [hm un]
  // sigma_hm_1_sqd_on_mu_hm_1_sqd: Ratio sigma_hm_1**2 / mu_hm_1**2    [-]
  // sigma_hm_2_sqd_on_mu_hm_2_sqd: Ratio sigma_hm_2**2 / mu_hm_2**2    [-]

  constexpr real_t mu_hm_min_coef = real_t(0.01); // Coef. for mu_hm_1_min and mu_hm_2_min

  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "calc_comp_mu_sigma_hm::cells_1", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          real_t Rmax{};           // Maximum possible value of ratio R                [-]
          real_t coef_A{};         // Coefficient A in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [-]
          real_t coef_B{};         // Coefficient B in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [hm un]
          real_t coef_C{};         // Coefficient C in A*mu_hm_1^2 + B*mu_hm_1 + C = 0 [hm un^2]
          real_t Bsqd_m_4AC{};     // Value B^2 - 4*A*C in quadratic eqn. for mu_hm_1  [hm un^2]
          real_t zeta_vrnce_rat{}; // Width parameter for sigma_hm_1^2 / mu_hm_1^2       [-]
          real_t mu_hm_1_min{};    // Minimum value of mu_hm_1 (precip. in both comps.) [hm un]
          real_t mu_hm_2_min{};    // Minimum value of mu_hm_2 (precip. in both comps.) [hm un]
          if (((hmm(i, k) >= hm_tol) && (precip_frac_1(i, k) >= precip_frac_tol(i)) &&
               (precip_frac_2(i, k) >= precip_frac_tol(i)))) {
            // Adjust the value of zeta based on the relationship of mu_thl_1 to
            // mu_thl_2.
            if ((mu_thl_1(i, k) <= mu_thl_2(i, k))) {
              if ((zeta_vrnce_rat_in(i) >= zero)) {
                zeta_vrnce_rat = zeta_vrnce_rat_in(i);
              } else {
                // zeta_vrnce_rat_in < 0
                zeta_vrnce_rat = ((one / (one + zeta_vrnce_rat_in(i))) - one);
              }
              // zeta_vrnce_rat_in >= 0
            } else {
              // mu_thl_1 > mu_thl_2
              if ((zeta_vrnce_rat_in(i) <= zero)) {
                zeta_vrnce_rat = zeta_vrnce_rat_in(i);
              } else {
                // zeta_vrnce_rat_in > 0
                zeta_vrnce_rat = ((one / (one + zeta_vrnce_rat_in(i))) - one);
              }
              // zeta_vrnce_rat_in <= 0
            }
            // mu_thl_1 <= mu_thl_2
            // Calculate the value of Rmax.
            // Rmax = ( f_p / ( a * f_p_1 * ( 1 + zeta ) + ( 1 - a ) * f_p_2 ) )
            // * ( <hm|_ip’^2> / <hm|_ip>^2 ).
            // The parameter zeta is written in the code as zeta_vrnce_rat.
            Rmax = ((precip_frac(i, k) /
                     (((mixt_frac(i, k) * precip_frac_1(i, k)) * (one + zeta_vrnce_rat)) +
                      ((one - mixt_frac(i, k)) * precip_frac_2(i, k)))) *
                    hmp2_ip_on_hmm2_ip);
            // Calculate the value of coefficient A.
            // A = a * f_p_1 * ( 1 + omicron * Rmax * ( 1 + zeta ) )
            // + a^2 * f_p_1^2 * ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ).
            coef_A = (((mixt_frac(i, k) * precip_frac_1(i, k)) *
                       (one + ((omicron(i) * Rmax) * (one + zeta_vrnce_rat)))) +
                      ((((mixt_frac(i, k) * mixt_frac(i, k)) *
                         (precip_frac_1(i, k) * precip_frac_1(i, k))) *
                        (one + (omicron(i) * Rmax))) /
                       ((one - mixt_frac(i, k)) * precip_frac_2(i, k))));
            // Calculate the value of coefficient B.
            // B = - 2 * <hm> * a * f_p_1 * ( 1 + omicron * Rmax )
            // / ( ( 1 - a ) * f_p_2 ).
            coef_B = (((((-(two)*hmm(i, k)) * mixt_frac(i, k)) * precip_frac_1(i, k)) *
                       (one + (omicron(i) * Rmax))) /
                      ((one - mixt_frac(i, k)) * precip_frac_2(i, k)));
            // Calculate the value of coefficient C.
            // C = - ( <hm’^2>
            // + ( 1 - ( 1 + omicron * Rmax ) / ( ( 1 - a ) * f_p_2 ) )
            // * <hm>^2 ).
            coef_C = -((hmp2(i, k) + ((one - ((one + (omicron(i) * Rmax)) /
                                              ((one - mixt_frac(i, k)) * precip_frac_2(i, k)))) *
                                      (hmm(i, k) * hmm(i, k)))));
            // Calculate value of B^2 - 4*A*C.
            Bsqd_m_4AC = ((coef_B * coef_B) - ((four * coef_A) * coef_C));
            // Mathematically, the value of B^2 - 4*A*C cannot be less than 0.
            // Numerically, this can happen when numerical round off error causes an
            // epsilon-sized negative value.  When this happens, reset the value of
            // B^2 - 4*A*C to 0.
            if ((Bsqd_m_4AC < zero)) {
              Bsqd_m_4AC = zero;
            }
            // Calculate the mean (in-precip.) of the hydrometeor in the 1st PDF
            // component.
            if ((mu_thl_1(i, k) <= mu_thl_2(i, k))) {
              mu_hm_1(i, k) = ((-(coef_B) + Kokkos::sqrt(Bsqd_m_4AC)) / (two * coef_A));
            } else {
              // mu_thl_1 > mu_thl_2
              mu_hm_1(i, k) = ((-(coef_B)-Kokkos::sqrt(Bsqd_m_4AC)) / (two * coef_A));
            }
            // mu_thl_1 <= mu_thl_2
            // Calculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
            // component.
            mu_hm_2(i, k) =
                ((hmm(i, k) - ((mixt_frac(i, k) * precip_frac_1(i, k)) * mu_hm_1(i, k))) /
                 ((one - mixt_frac(i, k)) * precip_frac_2(i, k)));
            // Calculate the value of the ratio R (which is sigma_hm_2^2 / mu_hm_2^2),
            // where R = omicron * Rmax.  The name of the variable used for R is
            // sigma_hm_2_sqd_on_mu_hm_2_sqd.
            sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) = (omicron(i) * Rmax);
            // Calculate minimum allowable values for mu_hm_1 and mu_hm_2.
            if (((hmm(i, k) / precip_frac(i, k)) > (hm_tol / precip_frac_1(i, k)))) {
              mu_hm_1_min =
                  (((hm_tol / precip_frac_1(i, k)) +
                    (mu_hm_min_coef *
                     ((hmm(i, k) / precip_frac(i, k)) - (hm_tol / precip_frac_1(i, k))))) <
                           ((hmm(i, k) - ((one - mixt_frac(i, k)) * hm_tol)) /
                            (mixt_frac(i, k) * precip_frac_1(i, k)))
                       ? ((hm_tol / precip_frac_1(i, k)) +
                          (mu_hm_min_coef *
                           ((hmm(i, k) / precip_frac(i, k)) - (hm_tol / precip_frac_1(i, k)))))
                       : ((hmm(i, k) - ((one - mixt_frac(i, k)) * hm_tol)) /
                          (mixt_frac(i, k) * precip_frac_1(i, k))));
            } else {
              // hmm / precip_frac <= hm_tol / precip_frac_1
              mu_hm_1_min = (hm_tol / precip_frac_1(i, k));
            }
            if (((hmm(i, k) / precip_frac(i, k)) > (hm_tol / precip_frac_2(i, k)))) {
              mu_hm_2_min =
                  (((hm_tol / precip_frac_2(i, k)) +
                    (mu_hm_min_coef *
                     ((hmm(i, k) / precip_frac(i, k)) - (hm_tol / precip_frac_2(i, k))))) <
                           ((hmm(i, k) - (mixt_frac(i, k) * hm_tol)) /
                            ((one - mixt_frac(i, k)) * precip_frac_2(i, k)))
                       ? ((hm_tol / precip_frac_2(i, k)) +
                          (mu_hm_min_coef *
                           ((hmm(i, k) / precip_frac(i, k)) - (hm_tol / precip_frac_2(i, k)))))
                       : ((hmm(i, k) - (mixt_frac(i, k) * hm_tol)) /
                          ((one - mixt_frac(i, k)) * precip_frac_2(i, k))));
            } else {
              // hmm / precip_frac <= hm_tol / precip_frac_2
              mu_hm_2_min = (hm_tol / precip_frac_2(i, k));
            }
            // Handle the "emergency" situation when the specified value of omicron is
            // too small for the value of <hm|_ip'^2> / <hm|_ip>^2, resulting in a
            // component mean that is too small (below tolerance value) or negative.
            if ((mu_hm_1(i, k) < mu_hm_1_min)) {
              // Set the value of mu_hm_1 to the threshold positive value.
              mu_hm_1(i, k) = mu_hm_1_min;
              // Recalculate the mean (in-precip.) of the hydrometeor in the 2nd PDF
              // component.
              mu_hm_2(i, k) =
                  ((hmm(i, k) - ((mixt_frac(i, k) * precip_frac_1(i, k)) * mu_hm_1(i, k))) /
                   ((one - mixt_frac(i, k)) * precip_frac_2(i, k)));
              // Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
              // scenario.
              // R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
              // - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
              // / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
              // + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
              // BFB-018: Fuse the left denominator product as in NVHPC Native CPU.
              sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) =
                  ((((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                     ((mixt_frac(i, k) * precip_frac_1(i, k)) * (mu_hm_1(i, k) * mu_hm_1(i, k)))) -
                    (((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                     (mu_hm_2(i, k) * mu_hm_2(i, k)))) /
                   Kokkos::fma((mixt_frac(i, k) * precip_frac_1(i, k)) * (one + zeta_vrnce_rat),
                               mu_hm_1(i, k) * mu_hm_1(i, k),
                               ((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                                   (mu_hm_2(i, k) * mu_hm_2(i, k))));
              // Mathematically, this ratio can never be less than 0.  In case
              // numerical round off error produces a negative value in extreme
              // cases, reset the value of R to 0.
              if ((sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) < zero)) {
                sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) = zero;
              }
            } else if ((mu_hm_2(i, k) < mu_hm_2_min)) {
              // Set the value of mu_hm_2 to the threshold positive value.
              mu_hm_2(i, k) = mu_hm_2_min;
              // Recalculate the mean (in-precip.) of the hydrometeor in the 1st PDF
              // component.
              mu_hm_1(i, k) =
                  ((hmm(i, k) - (((one - mixt_frac(i, k)) * precip_frac_2(i, k)) * mu_hm_2(i, k))) /
                   (mixt_frac(i, k) * precip_frac_1(i, k)));
              // Recalculate the value of R ( sigma_hm_2^2 / mu_hm_2^2 ) in this
              // scenario.
              // R = ( <hm'^2> + <hm>^2 - a * f_p_1 * mu_hm_1^2
              // - ( 1 - a ) * f_p_2 * mu_hm_2^2 )
              // / ( a * f_p_1 * ( 1 + zeta ) * mu_hm_1^2
              // + ( 1 - a ) * f_p_2 * mu_hm_2^2 ).
              // BFB-018: Fuse the left denominator product as in NVHPC Native CPU.
              sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) =
                  ((((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                     ((mixt_frac(i, k) * precip_frac_1(i, k)) * (mu_hm_1(i, k) * mu_hm_1(i, k)))) -
                    (((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                     (mu_hm_2(i, k) * mu_hm_2(i, k)))) /
                   Kokkos::fma((mixt_frac(i, k) * precip_frac_1(i, k)) * (one + zeta_vrnce_rat),
                               mu_hm_1(i, k) * mu_hm_1(i, k),
                               ((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                                   (mu_hm_2(i, k) * mu_hm_2(i, k))));
              // Mathematically, this ratio can never be less than 0.  In case
              // numerical round off error produces a negative value in extreme
              // cases, reset the value of R to 0.
              if ((sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) < zero)) {
                sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) = zero;
              }
            }
            // Calculate the standard deviation (in-precip.) of the hydrometeor in the
            // 1st PDF component.
            sigma_hm_1(i, k) =
                (Kokkos::sqrt((sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) * (one + zeta_vrnce_rat))) *
                 mu_hm_1(i, k));
            // Calculate the standard deviation (in-precip.) of the hydrometeor in the
            // 2nd PDF component.
            sigma_hm_2(i, k) = (Kokkos::sqrt(sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k)) * mu_hm_2(i, k));
            // Calculate the mean of the hydrometeor in the 1st PDF component.
            hm_1(i, k) = ((mu_hm_1(i, k) * precip_frac_1(i, k)) > hm_tol
                              ? (mu_hm_1(i, k) * precip_frac_1(i, k))
                              : hm_tol);
            // Calculate the mean of the hydrometeor in the 2nd PDF component.
            hm_2(i, k) = ((mu_hm_2(i, k) * precip_frac_2(i, k)) > hm_tol
                              ? (mu_hm_2(i, k) * precip_frac_2(i, k))
                              : hm_tol);
            // Calculate the ratio of sigma_hm_1^2 / mu_hm_1^2.
            sigma_hm_1_sqd_on_mu_hm_1_sqd(i, k) =
                ((sigma_hm_1(i, k) * sigma_hm_1(i, k)) / (mu_hm_1(i, k) * mu_hm_1(i, k)));
            // The value of R, sigma_hm_2_sqd_on_mu_hm_2_sqd, has already been
            // calculated.
          } else if (((hmm(i, k) >= hm_tol) && (precip_frac_1(i, k) >= precip_frac_tol(i)))) {
            // Precipitation is found in the 1st PDF component, but not in the 2nd
            // PDF component (precip_frac_2 = 0).
            mu_hm_1(i, k)    = (hmm(i, k) / (mixt_frac(i, k) * precip_frac_1(i, k)));
            mu_hm_2(i, k)    = zero;
            sigma_hm_1(i, k) = Kokkos::sqrt(
                ((((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                   ((mixt_frac(i, k) * precip_frac_1(i, k)) * (mu_hm_1(i, k) * mu_hm_1(i, k)))) /
                  (mixt_frac(i, k) * precip_frac_1(i, k))) > zero
                     ? (((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                         ((mixt_frac(i, k) * precip_frac_1(i, k)) *
                          (mu_hm_1(i, k) * mu_hm_1(i, k)))) /
                        (mixt_frac(i, k) * precip_frac_1(i, k)))
                     : zero));
            sigma_hm_2(i, k) = zero;
            hm_1(i, k)       = (mu_hm_1(i, k) * precip_frac_1(i, k));
            hm_2(i, k)       = zero;
            sigma_hm_1_sqd_on_mu_hm_1_sqd(i, k) =
                ((sigma_hm_1(i, k) * sigma_hm_1(i, k)) / (mu_hm_1(i, k) * mu_hm_1(i, k)));
            // The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
            sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) = zero;
          } else if (((hmm(i, k) >= hm_tol) && (precip_frac_2(i, k) >= precip_frac_tol(i)))) {
            // Precipitation is found in the 2nd PDF component, but not in the 1st
            // PDF component (precip_frac_1 = 0).
            mu_hm_1(i, k)    = zero;
            mu_hm_2(i, k)    = (hmm(i, k) / ((one - mixt_frac(i, k)) * precip_frac_2(i, k)));
            sigma_hm_1(i, k) = zero;
            sigma_hm_2(i, k) =
                Kokkos::sqrt(((((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                                (((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                                 (mu_hm_2(i, k) * mu_hm_2(i, k)))) /
                               ((one - mixt_frac(i, k)) * precip_frac_2(i, k))) > zero
                                  ? (((hmp2(i, k) + (hmm(i, k) * hmm(i, k))) -
                                      (((one - mixt_frac(i, k)) * precip_frac_2(i, k)) *
                                       (mu_hm_2(i, k) * mu_hm_2(i, k)))) /
                                     ((one - mixt_frac(i, k)) * precip_frac_2(i, k)))
                                  : zero));
            hm_1(i, k) = zero;
            hm_2(i, k) = (mu_hm_2(i, k) * precip_frac_2(i, k));
            // The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
            sigma_hm_1_sqd_on_mu_hm_1_sqd(i, k) = zero;
            sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) =
                ((sigma_hm_2(i, k) * sigma_hm_2(i, k)) / (mu_hm_2(i, k) * mu_hm_2(i, k)));
          } else {
            // hm < hm_tol or ( precip_frac_1 = 0 and precip_frac_2 = 0 ).
            // Precipitation is not found in either PDF component.
            mu_hm_1(i, k)    = zero;
            mu_hm_2(i, k)    = zero;
            sigma_hm_1(i, k) = zero;
            sigma_hm_2(i, k) = zero;
            hm_1(i, k)       = zero;
            hm_2(i, k)       = zero;
            // The ratio sigma_hm_1^2 / mu_hm_1^2 is undefined.
            sigma_hm_1_sqd_on_mu_hm_1_sqd(i, k) = zero;
            // The ratio sigma_hm_2^2 / mu_hm_2^2 is undefined.
            sigma_hm_2_sqd_on_mu_hm_2_sqd(i, k) = zero;
          }
          // hmm >= hm_tol and precip_frac_1 >= precip_frac_tol
          // and precip_frac_2 >= precip_frac_tol
        });
  }
  return;
}
// =============================================================================
namespace {
void component_corr_w_x(int nzt, int ngrdcol,                               // In
                        Real2D rc_1,                                        // In
                        Real2D rc_2,                                        // In
                        real_t corr_w_x_NN_cloud, real_t corr_w_x_NN_below, // In
                        int iiPDF_type,                                     // In
                        Real2D corr_w_x_1,                                  // Out
                        Real2D corr_w_x_2                                   // Out
) {
  // Description:
  // Calculates the correlation of w and x within the ith PDF component.
  // Here, x is a variable with a normally distributed individual marginal PDF,
  // such as chi or eta.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // rc_1: Mean cloud water mixing ratio (1st PDF comp.)  [kg/kg]
  // rc_2: Mean cloud water mixing ratio (2nd PDF comp.)  [kg/kg]
  // corr_w_x_NN_cloud: Corr. of w and x (ith PDF comp.); cloudy levs [-]
  // corr_w_x_NN_below: Corr. of w and x (ith PDF comp.); clear levs  [-]
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF.
  // corr_w_x_1: Correlation of w and x (1st PDF component)  [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF component)  [-]

  // Correlation of w and x in the ith PDF component.
  // The PDF variables chi and eta result from a transformation of the PDF
  // involving r_t and theta_l.  The correlation of w and x (whether x is chi
  // or eta) depends on the correlation of w and r_t, the correlation of w and
  // theta_l, as well as the variances of r_t and theta_l, and other factors.
  // The correlation of w and x is subject to change at every vertical level
  // and model time step, and is calculated as part of the CLUBB PDF
  // parameters.
  // The ADG1 PDF fixes the correlation of w and rt and the correlation of
  // w and theta_l to be 0, which means the correlation of w and chi and the
  // correlation of w and eta must also be 0.
  if ((((iiPDF_type == iiPDF_ADG1) || (iiPDF_type == iiPDF_ADG2) ||
        (iiPDF_type == iiPDF_new_hybrid)) &&
       l_follow_ADG1_PDF_standards)) {
    {
      const int slice_extent_0 = corr_w_x_1.extent_int(0);
      const int slice_extent_1 = corr_w_x_1.extent_int(1);
      Kokkos::parallel_for(
          "component_corr_w_x::array_1", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0            = ((cell % slice_extent_0));
            const int slice_1            = ((cell / slice_extent_0));
            corr_w_x_1(slice_0, slice_1) = zero;
          });
    }
    {
      const int slice_extent_0 = corr_w_x_2.extent_int(0);
      const int slice_extent_1 = corr_w_x_2.extent_int(1);
      Kokkos::parallel_for(
          "component_corr_w_x::array_2", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0            = ((cell % slice_extent_0));
            const int slice_1            = ((cell / slice_extent_0));
            corr_w_x_2(slice_0, slice_1) = zero;
          });
    }
  } else {
    // use prescribed paramter values
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "component_corr_w_x::cells_3", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            if ((rc_1(i, k) > rc_tol)) {
              corr_w_x_1(i, k) = corr_w_x_NN_cloud;
            } else {
              corr_w_x_1(i, k) = corr_w_x_NN_below;
            }
            if ((rc_2(i, k) > rc_tol)) {
              corr_w_x_2(i, k) = corr_w_x_NN_cloud;
            } else {
              corr_w_x_2(i, k) = corr_w_x_NN_below;
            }
          });
    }
  }
  // iiPDF_type
  return;
}
} // namespace
// =============================================================================
namespace {
void component_corr_chi_eta(int nzt, int ngrdcol,         // In
                            Real2D rc_1,                  // In
                            Real2D rc_2,                  // In
                            real_t corr_chi_eta_NN_cloud, // In
                            real_t corr_chi_eta_NN_below, // In
                            bool l_limit_corr_chi_eta,    // In
                            Real2D corr_chi_eta_1,        // Out
                            Real2D corr_chi_eta_2         // Out
) {
  // Description:
  // Calculates the correlation of chi (old s) and eta (old t) within the
  // ith PDF component.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // rc_1: Mean cloud water mix. rat. (1st PDF comp.) [kg/kg]&
  // rc_2: Mean cloud water mix. rat. (2nd PDF comp.) [kg/kg]
  // corr_chi_eta_NN_cloud: Corr. of chi & eta (ith PDF comp.); cloudy  [-]
  // corr_chi_eta_NN_below: Corr. of chi & eta (ith PDF comp.); clear   [-]
  // l_limit_corr_chi_eta: We must limit the correlation of chi and eta if we are to take the
  // Cholesky decomposition of the resulting correlation matrix. This is because a perfect
  // correlation of chi and eta was found to be unrealizable.
  // corr_chi_eta_1: Correlation of chi and eta (1st PDF component)     [-]
  // corr_chi_eta_2: Correlation of chi and eta (2nd PDF component)     [-]

  // Correlation of chi (old s) and eta (old t) in the ith PDF component.
  // The PDF variables chi and eta result from a transformation of the PDF
  // involving r_t and theta_l.  The correlation of chi and eta depends on the
  // correlation of r_t and theta_l, as well as the variances of r_t and
  // theta_l, and other factors.  The correlation of chi and eta is subject to
  // change at every vertical level and model time step, and is calculated as
  // part of the CLUBB PDF parameters.
  // WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
  // code is necessary because SILHS is lazy and wussy, and only
  // wants to declare correlation arrays at the start of the model
  // run, rather than updating them throughout the model run.
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "component_corr_chi_eta::cells_1", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((rc_1(i, k) > rc_tol)) {
            corr_chi_eta_1(i, k) = corr_chi_eta_NN_cloud;
          } else {
            corr_chi_eta_1(i, k) = corr_chi_eta_NN_below;
          }
        });
  }
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "component_corr_chi_eta::cells_2", ((extent_0 * extent_1)), KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((rc_2(i, k) > rc_tol)) {
            corr_chi_eta_2(i, k) = corr_chi_eta_NN_cloud;
          } else {
            corr_chi_eta_2(i, k) = corr_chi_eta_NN_below;
          }
        });
  }
  // We cannot have a perfect correlation of chi (old s) and eta (old t) if we
  // plan to decompose this matrix and we don't want the Cholesky_factor code
  // to throw a fit.
  if (l_limit_corr_chi_eta) {
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "component_corr_chi_eta::cells_3", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            corr_chi_eta_1(i, k) =
                ((corr_chi_eta_1(i, k) < max_mag_correlation
                      ? corr_chi_eta_1(i, k)
                      : max_mag_correlation) > -(max_mag_correlation)
                     ? (corr_chi_eta_1(i, k) < max_mag_correlation ? corr_chi_eta_1(i, k)
                                                                   : max_mag_correlation)
                     : -(max_mag_correlation));
          });
    }
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "component_corr_chi_eta::cells_4", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            corr_chi_eta_2(i, k) =
                ((corr_chi_eta_2(i, k) < max_mag_correlation
                      ? corr_chi_eta_2(i, k)
                      : max_mag_correlation) > -(max_mag_correlation)
                     ? (corr_chi_eta_2(i, k) < max_mag_correlation ? corr_chi_eta_2(i, k)
                                                                   : max_mag_correlation)
                     : -(max_mag_correlation));
          });
    }
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void component_corr_w_hm_n_ip(int nzt, int ngrdcol,                 // In
                              Real2D corr_w_hm_1_n_in, Real2D rc_1, // In
                              Real2D corr_w_hm_2_n_in, Real2D rc_2, // In
                              real_t corr_w_hm_n_NL_cloud,          // In
                              real_t corr_w_hm_n_NL_below,          // In
                              bool l_calc_w_corr,                   // In
                              Real2D corr_w_hm_1_n,                 // Out
                              Real2D corr_w_hm_2_n                  // Out
) {
  // Description:
  // Calculates the in-precip correlation of w and the natural logarithm of a
  // hydrometeor species within the ith PDF component.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // corr_w_hm_1_n_in: Correlation of w and ln hm (1st PDF comp.) ip    [-]
  // rc_1: Mean cloud water mix. ratio (1st PDF comp.)  [kg/kg]
  // corr_w_hm_2_n_in: Correlation of w and ln hm (2nd PDF comp.) ip    [-]
  // rc_2: Mean cloud water mix. ratio (2nd PDF comp.)  [kg/kg]
  // corr_w_hm_n_NL_cloud: Corr. of w & ln hm (ith PDF comp.) ip; cloud [-]
  // corr_w_hm_n_NL_below: Corr. of w & ln hm (ith PDF comp.) ip; clear [-]
  // l_calc_w_corr: Calculate the correlations between w and the hydrometeors
  // corr_w_hm_1_n: Correlation of w and ln hm (1st PDF component) ip   [-]
  // corr_w_hm_2_n: Correlation of w and ln hm (2nd PDF component) ip   [-]

  // Correlation (in-precip) of w and the natural logarithm of the hydrometeor
  // in the ith PDF component.
  if (l_calc_w_corr) {
    {
      const int slice_extent_0 = corr_w_hm_1_n.extent_int(0);
      const int slice_extent_1 = corr_w_hm_1_n.extent_int(1);
      Kokkos::parallel_for(
          "component_corr_w_hm_n_ip::array_1", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = ((cell % slice_extent_0));
            const int slice_1               = ((cell / slice_extent_0));
            corr_w_hm_1_n(slice_0, slice_1) = corr_w_hm_1_n_in(slice_0, slice_1);
          });
    }
    {
      const int slice_extent_0 = corr_w_hm_2_n.extent_int(0);
      const int slice_extent_1 = corr_w_hm_2_n.extent_int(1);
      Kokkos::parallel_for(
          "component_corr_w_hm_n_ip::array_2", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = ((cell % slice_extent_0));
            const int slice_1               = ((cell / slice_extent_0));
            corr_w_hm_2_n(slice_0, slice_1) = corr_w_hm_2_n_in(slice_0, slice_1);
          });
    }
  } else {
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "component_corr_w_hm_n_ip::cells_3", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            if ((rc_1(i, k) > rc_tol)) {
              corr_w_hm_1_n(i, k) = corr_w_hm_n_NL_cloud;
            } else {
              corr_w_hm_1_n(i, k) = corr_w_hm_n_NL_below;
            }
            if ((rc_2(i, k) > rc_tol)) {
              corr_w_hm_2_n(i, k) = corr_w_hm_n_NL_cloud;
            } else {
              corr_w_hm_2_n(i, k) = corr_w_hm_n_NL_below;
            }
          });
    }
  }
  // l_calc_w_corr
  return;
}
} // namespace
// =============================================================================
namespace {
void component_corr_x_hm_n_ip(int nzt, int ngrdcol,        // In
                              Real2D rc_1,                 // In
                              Real2D rc_2,                 // In
                              real_t corr_x_hm_n_NL_cloud, // In
                              real_t corr_x_hm_n_NL_below, // In
                              Real2D corr_x_hm_1_n,        // Out
                              Real2D corr_x_hm_2_n         // Out
) {
  // Description:
  // Calculates the in-precip correlation of x and a hydrometeor species
  // within the ith PDF component.  Here, x is a variable with a normally
  // distributed individual marginal PDF, such as chi or eta.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // rc_1: Mean cloud water mixing ratio (1st PDF comp.)   [kg/kg]
  // rc_2: Mean cloud water mixing ratio (2nd PDF comp.)   [kg/kg]
  // corr_x_hm_n_NL_cloud: Corr. of x and ln hm (ith PDF comp.) ip     [-]
  // corr_x_hm_n_NL_below: Corr. of x and ln hm (ith PDF comp.) ip     [-]
  // corr_x_hm_1_n: Correlation of x and ln hm (1st PDF component) ip  [-]
  // corr_x_hm_2_n: Correlation of x and ln hm (2nd PDF component) ip  [-]

  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "component_corr_x_hm_n_ip::cells_1", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((rc_1(i, k) > rc_tol)) {
            corr_x_hm_1_n(i, k) = corr_x_hm_n_NL_cloud;
          } else {
            corr_x_hm_1_n(i, k) = corr_x_hm_n_NL_below;
          }
          if ((rc_2(i, k) > rc_tol)) {
            corr_x_hm_2_n(i, k) = corr_x_hm_n_NL_cloud;
          } else {
            corr_x_hm_2_n(i, k) = corr_x_hm_n_NL_below;
          }
        });
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void component_corr_hmx_hmy_n_ip(int nzt, int ngrdcol,           // In
                                 Real2D rc_1,                    // In
                                 Real2D rc_2,                    // In
                                 real_t corr_hmx_hmy_n_LL_cloud, // In
                                 real_t corr_hmx_hmy_n_LL_below, // In
                                 Real2D corr_hmx_hmy_1_n,        // Out
                                 Real2D corr_hmx_hmy_2_n         // Out
) {
  // Description:
  // Calculates the in-precip correlation of the natural logarithms of
  // hydrometeor x and hydrometeor y within the ith PDF component.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // rc_1: Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
  // rc_2: Mean cloud water mixing ratio (ith PDF comp.) [kg/kg]
  // corr_hmx_hmy_n_LL_cloud: Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]
  // corr_hmx_hmy_n_LL_below: Corr.: ln hmx & ln hmy (ith PDF comp.) ip [-]
  // corr_hmx_hmy_1_n: Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]
  // corr_hmx_hmy_2_n: Corr. of ln hmx & ln hmy (ith PDF comp.) ip      [-]

  // Correlation (in-precip) of the natural logarithms of hydrometeor x and
  // hydrometeor y in the ith PDF component.
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "component_corr_hmx_hmy_n_ip::cells_1", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((rc_1(i, k) > rc_tol)) {
            corr_hmx_hmy_1_n(i, k) = corr_hmx_hmy_n_LL_cloud;
          } else {
            corr_hmx_hmy_1_n(i, k) = corr_hmx_hmy_n_LL_below;
          }
          if ((rc_2(i, k) > rc_tol)) {
            corr_hmx_hmy_2_n(i, k) = corr_hmx_hmy_n_LL_cloud;
          } else {
            corr_hmx_hmy_2_n(i, k) = corr_hmx_hmy_n_LL_below;
          }
        });
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void component_corr_eta_hm_n_ip(int nzt, int ngrdcol,   // In
                                Real2D corr_chi_eta_1,  // In
                                Real2D corr_chi_hm_n_1, // In
                                Real2D corr_chi_eta_2,  // In
                                Real2D corr_chi_hm_n_2, // In
                                Real2D corr_eta_hm_n_1, // Out
                                Real2D corr_eta_hm_n_2  // Out
) {
  // Description:
  // Estimates the correlation of eta and the natural logarithm of a
  // hydrometeor species using the correlation of chi and eta and the
  // correlation of chi and the natural logarithm of the hydrometeor.  This
  // facilities the Cholesky decomposability of the correlation array that will
  // inevitably be decomposed for SILHS purposes. Without this estimation, we
  // have found that the resulting correlation matrix cannot be decomposed.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // corr_chi_eta_1: Component correlation of chi and eta              [-]
  // corr_chi_hm_n_1: Component correlation of chi and ln hm            [-]
  // corr_chi_eta_2: Component correlation of chi and eta              [-]
  // corr_chi_hm_n_2: Component correlation of chi and ln hm            [-]
  // corr_eta_hm_n_1: Component correlation of eta and ln hm            [-]
  // corr_eta_hm_n_2: Component correlation of eta and ln hm            [-]

  {
    const int slice_extent_0 = corr_eta_hm_n_1.extent_int(0);
    const int slice_extent_1 = corr_eta_hm_n_1.extent_int(1);
    Kokkos::parallel_for(
        "component_corr_eta_hm_n_ip::array_1", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_eta_hm_n_1(slice_0, slice_1) =
              (corr_chi_eta_1(slice_0, slice_1) * corr_chi_hm_n_1(slice_0, slice_1));
        });
  }
  {
    const int slice_extent_0 = corr_eta_hm_n_2.extent_int(0);
    const int slice_extent_1 = corr_eta_hm_n_2.extent_int(1);
    Kokkos::parallel_for(
        "component_corr_eta_hm_n_ip::array_2", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_eta_hm_n_2(slice_0, slice_1) =
              (corr_chi_eta_2(slice_0, slice_1) * corr_chi_hm_n_2(slice_0, slice_1));
        });
  }
  return;
}
} // namespace
// =============================================================================
void norm_transform_mean_stdev(int nzt, int ngrdcol,                                  // In
                               int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                               Real3D hm_1, Real3D hm_2,                              // In
                               Real2D Ncnm, int pdf_dim,                              // In
                               Real3D mu_x_1, Real3D mu_x_2,                          // In
                               Real3D sigma_x_1, Real3D sigma_x_2,                    // In
                               Real3D sigma2_on_mu2_ip_1,                             // In
                               Real3D sigma2_on_mu2_ip_2,                             // In
                               bool l_const_Nc_in_cloud,                              // In
                               Real3D mu_x_1_n, Real3D mu_x_2_n,                      // Out
                               Real3D sigma_x_1_n, Real3D sigma_x_2_n                 // Out
) {
  // Description:
  // Transforms the means and the standard deviations of PDF variables that
  // have assumed lognormal distributions -- which are precipitating
  // hydrometeors (in precipitation) and N_cn -- to normal space for each PDF
  // component.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // hm_1: Mean of a precip. hydrometeor (1st PDF component)   [units vary]
  // hm_2: Mean of a precip. hydrometeor (2nd PDF component)   [units vary]
  // Ncnm: Mean cloud nuclei concentration, < N_cn >               [num/kg]
  // pdf_dim: Number of variables in CLUBB's PDF
  // mu_x_1: Mean array of PDF vars. (1st PDF component)    [units vary]
  // mu_x_2: Mean array of PDF vars. (2nd PDF component)    [units vary]
  // sigma_x_1: Standard deviation array of PDF vars (comp. 1) [units vary]
  // sigma_x_2: Standard deviation array of PDF vars (comp. 2) [units vary]
  // sigma2_on_mu2_ip_1: Prescribed ratio array: sigma_hm_1^2/mu_hm_1^2 [-]
  // sigma2_on_mu2_ip_2: Prescribed ratio array: sigma_hm_2^2/mu_hm_2^2 [-]
  // l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K)
  // mu_x_1_n: Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // mu_x_2_n: Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

  int ivar{};
  int hm_idx{};
  int k{};
  int iiPDF_w{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi    = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta    = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w      = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn    = hm_metadata.iiPDF_Ncn;
  const auto hm_metadata_hydromet_tol = hm_metadata.hydromet_tol;

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // The means and standard deviations in each PDF component of w, chi (old s),
  // and eta (old t) do not need to be transformed to normal space, since w,
  // chi, and eta already follow assumed normal distributions in each PDF
  // component.  The normal space means and standard deviations are the same as
  // the actual means and standard deviations.
  for (k = 0; k < nzt; ++k) {
    {
      const int slice_extent_0 = mu_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_1", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = cell;
            mu_x_1_n(slice_0, k, iiPDF_chi) = mu_x_1(slice_0, k, iiPDF_chi);
          });
    }
    {
      const int slice_extent_0 = mu_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_2", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = cell;
            mu_x_1_n(slice_0, k, iiPDF_eta) = mu_x_1(slice_0, k, iiPDF_eta);
          });
    }
    {
      const int slice_extent_0 = mu_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_3", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0             = cell;
            mu_x_1_n(slice_0, k, iiPDF_w) = mu_x_1(slice_0, k, iiPDF_w);
          });
    }
  }
  for (k = 0; k < nzt; ++k) {
    {
      const int slice_extent_0 = mu_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_4", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = cell;
            mu_x_2_n(slice_0, k, iiPDF_chi) = mu_x_2(slice_0, k, iiPDF_chi);
          });
    }
    {
      const int slice_extent_0 = mu_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_5", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0               = cell;
            mu_x_2_n(slice_0, k, iiPDF_eta) = mu_x_2(slice_0, k, iiPDF_eta);
          });
    }
    {
      const int slice_extent_0 = mu_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_6", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0             = cell;
            mu_x_2_n(slice_0, k, iiPDF_w) = mu_x_2(slice_0, k, iiPDF_w);
          });
    }
  }
  for (k = 0; k < nzt; ++k) {
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_7", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                  = cell;
            sigma_x_1_n(slice_0, k, iiPDF_chi) = sigma_x_1(slice_0, k, iiPDF_chi);
          });
    }
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_8", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                  = cell;
            sigma_x_1_n(slice_0, k, iiPDF_eta) = sigma_x_1(slice_0, k, iiPDF_eta);
          });
    }
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_9", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                = cell;
            sigma_x_1_n(slice_0, k, iiPDF_w) = sigma_x_1(slice_0, k, iiPDF_w);
          });
    }
  }
  for (k = 0; k < nzt; ++k) {
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_10", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                  = cell;
            sigma_x_2_n(slice_0, k, iiPDF_chi) = sigma_x_2(slice_0, k, iiPDF_chi);
          });
    }
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_11", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                  = cell;
            sigma_x_2_n(slice_0, k, iiPDF_eta) = sigma_x_2(slice_0, k, iiPDF_eta);
          });
    }
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_12", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0                = cell;
            sigma_x_2_n(slice_0, k, iiPDF_w) = sigma_x_2(slice_0, k, iiPDF_w);
          });
    }
  }
  // !! Transform the mean and standard deviation to normal space in each PDF
  // !! component for variables that have an assumed lognormal distribution,
  // !! given the mean and standard deviation in each PDF component for those
  // !! variables.  A precipitating hydrometeor has an assumed lognormal
  // !! distribution in precipitation in each PDF component.  Simplified cloud
  // !! nuclei concentration, N_cn, has an assumed lognormal distribution in
  // !! each PDF component, and furthermore, mu_Ncn_1 = mu_Ncn_2 and
  // !! sigma_Ncn_1 = sigma_Ncn_2, so N_cn has an assumed single lognormal
  // !! distribution over the entire domain.
  // Normal space mean of simplified cloud nuclei concentration, N_cn,
  // in PDF component 1.
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "norm_transform_mean_stdev::cells_13", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((Ncnm(i, k) >= Ncn_tol)) {
            mu_x_1_n(i, k, iiPDF_Ncn) =
                mean_L2N(mu_x_1(i, k, iiPDF_Ncn), sigma2_on_mu2_ip_1(i, k, iiPDF_Ncn));
          } else {
            // Mean simplified cloud nuclei concentration in PDF component 1 is less
            // than the tolerance amount.  It is considered to have a value of 0.
            // There are not any cloud nuclei or cloud at this grid level.  The value
            // of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
            // assigning it a value.
            mu_x_1_n(i, k, iiPDF_Ncn) = -(Kokkos::finite_max_v<real_t>);
          }
        });
  }
  // Normal space mean of simplified cloud nuclei concentration, N_cn,
  // in PDF component 2.
  {
    const int extent_0 = nzt;
    const int extent_1 = ngrdcol;
    Kokkos::parallel_for(
        "norm_transform_mean_stdev::cells_14", ((extent_0 * extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % extent_1));
          const int k = ((cell / extent_1));
          if ((Ncnm(i, k) >= Ncn_tol)) {
            mu_x_2_n(i, k, iiPDF_Ncn) =
                mean_L2N(mu_x_2(i, k, iiPDF_Ncn), sigma2_on_mu2_ip_1(i, k, iiPDF_Ncn));
          } else {
            // Mean simplified cloud nuclei concentration in PDF component 1 is less
            // than the tolerance amount.  It is considered to have a value of 0.
            // There are not any cloud nuclei or cloud at this grid level.  The value
            // of mu_Ncn_1_n should be -inf.  It will be set to -huge for purposes of
            // assigning it a value.
            mu_x_2_n(i, k, iiPDF_Ncn) = -(Kokkos::finite_max_v<real_t>);
          }
        });
  }
  // Normal space standard deviation of simplified cloud nuclei concentration,
  // N_cn, in PDF components 1 and 2.
  if (l_const_Nc_in_cloud) {
    // Ncn does not vary in the grid box.
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      const int slice_extent_1 = sigma_x_1_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_15", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                        = ((cell % slice_extent_0));
            const int slice_1                        = ((cell / slice_extent_0));
            sigma_x_1_n(slice_0, slice_1, iiPDF_Ncn) = zero;
          });
    }
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      const int slice_extent_1 = sigma_x_2_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_16", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                        = ((cell % slice_extent_0));
            const int slice_1                        = ((cell / slice_extent_0));
            sigma_x_2_n(slice_0, slice_1, iiPDF_Ncn) = zero;
          });
    }
  } else {
    // Ncn (perhaps) varies in the grid box.
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      const int slice_extent_1 = sigma_x_1_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_17", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_1_n(slice_0, slice_1, iiPDF_Ncn) =
                stdev_L2N(sigma2_on_mu2_ip_1(slice_0, slice_1, iiPDF_Ncn));
          });
    }
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      const int slice_extent_1 = sigma_x_2_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_18", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_2_n(slice_0, slice_1, iiPDF_Ncn) =
                stdev_L2N(sigma2_on_mu2_ip_2(slice_0, slice_1, iiPDF_Ncn));
          });
    }
  }
  // Normal space precipitating hydrometeor means and standard deviations.
  for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
    hm_idx = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
    // Normal space mean of a precipitating hydrometeor, hm, in PDF
    // component 1.
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::cells_19", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            if ((hm_1(i, k, hm_idx) >= hm_metadata_hydromet_tol(hm_idx))) {
              mu_x_1_n(i, k, ivar) = mean_L2N(mu_x_1(i, k, ivar), sigma2_on_mu2_ip_1(i, k, ivar));
            } else {
              // The mean of a precipitating hydrometeor in PDF component 1 is less
              // than its tolerance amount.  It is considered to have a value of 0.
              // There is not any of this precipitating hydrometeor in the 1st PDF
              // component at this grid level.  The in-precip mean of this
              // precipitating hydrometeor (1st PDF component) is also 0.  The value
              // of mu_hm_1_n should be -inf.  It will be set to -huge for purposes
              // of assigning it a value.
              mu_x_1_n(i, k, ivar) = -(Kokkos::finite_max_v<real_t>);
            }
          });
    }
    // Normal space standard deviation of a precipitating hydrometeor, hm, in
    // PDF component 1.
    {
      const int slice_extent_0 = sigma_x_1_n.extent_int(0);
      const int slice_extent_1 = sigma_x_1_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_20", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_1_n(slice_0, slice_1, ivar) =
                stdev_L2N(sigma2_on_mu2_ip_1(slice_0, slice_1, ivar));
          });
    }
    // Normal space mean of a precipitating hydrometeor, hm, in PDF
    // component 2.
    {
      const int extent_0 = nzt;
      const int extent_1 = ngrdcol;
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::cells_21", ((extent_0 * extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % extent_1));
            const int k = ((cell / extent_1));
            if ((hm_2(i, k, hm_idx) >= hm_metadata_hydromet_tol(hm_idx))) {
              mu_x_2_n(i, k, ivar) = mean_L2N(mu_x_2(i, k, ivar), sigma2_on_mu2_ip_2(i, k, ivar));
            } else {
              // The mean of a precipitating hydrometeor in PDF component 2 is less
              // than its tolerance amount.  It is considered to have a value of 0.
              // There is not any of this precipitating hydrometeor in the 2nd PDF
              // component at this grid level.  The in-precip mean of this
              // precipitating hydrometeor (2nd PDF component) is also 0.  The value
              // of mu_hm_2_n should be -inf.  It will be set to -huge for purposes
              // of assigning it a value.
              mu_x_2_n(i, k, ivar) = -(Kokkos::finite_max_v<real_t>);
            }
          });
    }
    // Normal space standard deviation of a precipitating hydrometeor, hm, in
    // PDF component 2.
    {
      const int slice_extent_0 = sigma_x_2_n.extent_int(0);
      const int slice_extent_1 = sigma_x_2_n.extent_int(1);
      Kokkos::parallel_for(
          "norm_transform_mean_stdev::array_22", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            sigma_x_2_n(slice_0, slice_1, ivar) =
                stdev_L2N(sigma2_on_mu2_ip_2(slice_0, slice_1, ivar));
          });
    }
  }
  // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
  return;
}
// =============================================================================
void denorm_transform_corr(int nzt, int ngrdcol, int pdf_dim,
                           const hm_metadata_type& hm_metadata,                  // In
                           Real3D sigma_x_1_n, Real3D sigma_x_2_n,               // In
                           Real3D sigma2_on_mu2_ip_1, Real3D sigma2_on_mu2_ip_2, // In
                           Real4D corr_array_1_n,                                // In
                           Real4D corr_array_2_n,                                // In
                           Real4D corr_array_1, Real4D corr_array_2              // Out
) {
  // Description:
  // Calculates the true or "real-space" correlations between PDF variables,
  // where at least one of the variables that is part of a correlation has an
  // assumed lognormal distribution -- which are the precipitating hydrometeors
  // (in precipitation) and N_cn.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical levels
  // ngrdcol: Number of grid columns
  // pdf_dim: Number of variables
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // sigma2_on_mu2_ip_1: Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
  // sigma2_on_mu2_ip_2: Ratio array sigma_hm_2^2/mu_hm_2^2             [-]
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // corr_array_1: Correlation array of PDF vars. (comp. 1)             [-]
  // corr_array_2: Correlation array of PDF vars. (comp. 2)             [-]

  int ivar{};
  int jvar{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_w{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w   = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn = hm_metadata.iiPDF_Ncn;

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // Native writes the diagonal and lower triangle only. Upper entries are
  // undefined; do not import the JAX symmetrization into this transform.
  // Initialize diagonal elements to one
  for (ivar = 0; ivar < pdf_dim; ++ivar) {
    {
      const int slice_extent_0 = corr_array_1.extent_int(0);
      const int slice_extent_1 = corr_array_1.extent_int(1);
      Kokkos::parallel_for(
          "denorm_transform_corr::array_1", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                          = ((cell % slice_extent_0));
            const int slice_1                          = ((cell / slice_extent_0));
            corr_array_1(slice_0, slice_1, ivar, ivar) = one;
          });
    }
    {
      const int slice_extent_0 = corr_array_2.extent_int(0);
      const int slice_extent_1 = corr_array_2.extent_int(1);
      Kokkos::parallel_for(
          "denorm_transform_corr::array_2", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                          = ((cell % slice_extent_0));
            const int slice_1                          = ((cell / slice_extent_0));
            corr_array_2(slice_0, slice_1, ivar, ivar) = one;
          });
    }
  }
  // The correlations in each PDF component between two of w, chi (old s), and
  // eta (old t) do not need to be transformed to standard space, since w, chi,
  // and eta follow assumed normal distributions in each PDF component.  The
  // normal space correlations between any two of these variables are the same
  // as the actual correlations.
  {
    const int slice_extent_0 = corr_array_1.extent_int(0);
    const int slice_extent_1 = corr_array_1.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_3", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_1(slice_0, slice_1, iiPDF_eta, iiPDF_chi) =
              corr_array_1_n(slice_0, slice_1, iiPDF_eta, iiPDF_chi);
        });
  }
  {
    const int slice_extent_0 = corr_array_1.extent_int(0);
    const int slice_extent_1 = corr_array_1.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_4", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_1(slice_0, slice_1, iiPDF_w, iiPDF_chi) =
              corr_array_1_n(slice_0, slice_1, iiPDF_w, iiPDF_chi);
        });
  }
  {
    const int slice_extent_0 = corr_array_1.extent_int(0);
    const int slice_extent_1 = corr_array_1.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_5", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_1(slice_0, slice_1, iiPDF_w, iiPDF_eta) =
              corr_array_1_n(slice_0, slice_1, iiPDF_w, iiPDF_eta);
        });
  }
  {
    const int slice_extent_0 = corr_array_2.extent_int(0);
    const int slice_extent_1 = corr_array_2.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_6", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_2(slice_0, slice_1, iiPDF_eta, iiPDF_chi) =
              corr_array_2_n(slice_0, slice_1, iiPDF_eta, iiPDF_chi);
        });
  }
  {
    const int slice_extent_0 = corr_array_2.extent_int(0);
    const int slice_extent_1 = corr_array_2.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_7", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_2(slice_0, slice_1, iiPDF_w, iiPDF_chi) =
              corr_array_2_n(slice_0, slice_1, iiPDF_w, iiPDF_chi);
        });
  }
  {
    const int slice_extent_0 = corr_array_2.extent_int(0);
    const int slice_extent_1 = corr_array_2.extent_int(1);
    Kokkos::parallel_for(
        "denorm_transform_corr::array_8", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          corr_array_2(slice_0, slice_1, iiPDF_w, iiPDF_eta) =
              corr_array_2_n(slice_0, slice_1, iiPDF_w, iiPDF_eta);
        });
  }
  // !! Calculate the true correlation of variables that have an assumed normal
  // !! distribution and variables that have an assumed lognormal distribution
  // !! for the ith PDF component, given their normal space correlation and the
  // !! normal space standard deviation of the variable with the assumed
  // !! lognormal distribution.
  // Transform the correlations between chi/eta/w and N_cn to standard space.
  // Transform the correlation of chi (old s) and N_cn to standard space in PDF component 1.
  corr_NN2NL(nzt, ngrdcol,                                                                    // In
             Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi), // In
             Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),               // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),        // In
             Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi)    // Out
  );
  // Transform the correlation of eta (old t) and N_cn to standard space in PDF component 1.
  corr_NN2NL(nzt, ngrdcol,                                                                    // In
             Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta), // In
             Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),               // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),        // In
             Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta)    // Out
  );
  // Transform the correlation of w and N_cn to standard space in PDF component 1.
  corr_NN2NL(nzt, ngrdcol,                                                                  // In
             Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
             Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),             // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),      // In
             Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w)    // Out
  );
  // Transform the correlation of chi (old s) and N_cn to standard space in PDF component 2.
  corr_NN2NL(nzt, ngrdcol,                                                                    // In
             Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi), // In
             Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),               // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),        // In
             Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi)    // Out
  );
  // Transform the correlation of eta (old t) and N_cn to standard space in PDF component 2.
  corr_NN2NL(nzt, ngrdcol,                                                                    // In
             Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta), // In
             Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),               // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),        // In
             Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta)    // Out
  );
  // Transform the correlation of w and N_cn to standard space in PDF component 2.
  corr_NN2NL(nzt, ngrdcol,                                                                  // In
             Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
             Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),             // In
             Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),      // In
             Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w)    // Out
  );
  // Transform the correlations (in-precip) between chi/eta/w and the
  // precipitating hydrometeors to standard space.
  for (ivar = iiPDF_chi; ivar < iiPDF_w + 1; ++ivar) {
    for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
      // Transform the correlation (in-precip) between w, chi, or eta and a
      // precipitating hydrometeor, hm, to standard space in PDF component 1.
      corr_NN2NL(nzt, ngrdcol,                                                          // In
                 Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                 Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, jvar),          // In
                 Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, jvar),   // In
                 Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, jvar, ivar)    // Out
      );
      // Transform the correlation (in-precip) between w, chi, or eta and a
      // precipitating hydrometeor, hm, to standard space in PDF component 2.
      corr_NN2NL(nzt, ngrdcol,                                                          // In
                 Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                 Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, jvar),          // In
                 Kokkos::subview(sigma2_on_mu2_ip_2, Kokkos::ALL, Kokkos::ALL, jvar),   // In
                 Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, jvar, ivar)    // Out
      );
    }
    // Native loop: jvar = iiPDF_Ncn+1, pdf_dim
  }
  // Native loop: ivar = iiPDF_chi, iiPDF_w
  // !! Calculate the true correlation of two variables that both have an
  // !! assumed lognormal distribution for the ith PDF component, given their
  // !! normal space correlation and both of their normal space standard
  // !! deviations.
  // Transform the correlations (in-precip) between N_cn and the precipitating
  // hydrometeors to standard space.
  for (jvar = ((iiPDF_Ncn + 1)); jvar < pdf_dim; ++jvar) {
    // Transform the correlation (in-precip) between N_cn and a precipitating
    // hydrometeor, hm, to standard space in PDF component 1.
    corr_NN2LL(nzt, ngrdcol,                                                               // In
               Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn), // In
               Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),          // In
               Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, jvar),               // In
               Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),   // In
               Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, jvar),        // In
               Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn)    // Out
    );
    // Transform the correlation (in-precip) between N_cn and a precipitating
    // hydrometeor, hm, to standard space in PDF component 2.
    corr_NN2LL(nzt, ngrdcol,                                                               // In
               Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn), // In
               Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),          // In
               Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, jvar),               // In
               Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn),   // In
               Kokkos::subview(sigma2_on_mu2_ip_2, Kokkos::ALL, Kokkos::ALL, jvar),        // In
               Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, jvar, iiPDF_Ncn)    // Out
    );
  }
  // Native loop: jvar = ivar+1, pdf_dim
  // Transform the correlations (in-precip) between two precipitating
  // hydrometeors to standard space.
  for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim - 1; ++ivar) {
    for (jvar = ((ivar + 1)); jvar < pdf_dim; ++jvar) {
      // Transform the correlation (in-precip) between two precipitating
      // hydrometeors (for example, r_r and N_r) to standard space in PDF
      // component 1.
      corr_NN2LL(nzt, ngrdcol,                                                          // In
                 Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                 Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, ivar),          // In
                 Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, jvar),          // In
                 Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, ivar),   // In
                 Kokkos::subview(sigma2_on_mu2_ip_1, Kokkos::ALL, Kokkos::ALL, jvar),   // In
                 Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, jvar, ivar)    // Out
      );
      // Transform the correlation (in-precip) between two precipitating
      // hydrometeors (for example, r_r and N_r) to standard space in PDF
      // component 2.
      corr_NN2LL(nzt, ngrdcol,                                                          // In
                 Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                 Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, ivar),          // In
                 Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, jvar),          // In
                 Kokkos::subview(sigma2_on_mu2_ip_2, Kokkos::ALL, Kokkos::ALL, ivar),   // In
                 Kokkos::subview(sigma2_on_mu2_ip_2, Kokkos::ALL, Kokkos::ALL, jvar),   // In
                 Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, jvar, ivar)    // Out
      );
    }
    // Native loop: jvar = ivar+1, pdf_dim
  }
  // Native loop: ivar = iiPDF_Ncn+1, pdf_dim-1
  return;
}
// =============================================================================
namespace {
KOKKOS_INLINE_FUNCTION void calc_corr_w_hm_n(real_t wm, real_t wphydrometp,               // In
                                             real_t mu_w_1, real_t mu_w_2,                // In
                                             real_t mu_hm_1, real_t mu_hm_2,              // In
                                             real_t sigma_w_1, real_t sigma_w_2,          // In
                                             real_t sigma_hm_1, real_t sigma_hm_2,        // In
                                             real_t sigma_hm_1_n, real_t sigma_hm_2_n,    // In
                                             real_t mixt_frac, real_t precip_frac_1,      // In
                                             real_t precip_frac_2, real_t hm_tol,         // In
                                             real_t& corr_w_hm_1_n, real_t& corr_w_hm_2_n // Out
) {
  // Description:
  // Calculates the PDF component correlation (in-precip) between vertical
  // velocity, w, and the natural logarithm of a hydrometeor, ln hm.  The
  // overall covariance of w and hm, <w'hm'> can be written in terms of the PDF
  // parameters.  When both w and hm vary in both PDF components, the equation
  // is written as:
  // <w'hm'> = mixt_frac * precip_frac_1
  //           * ( mu_w_1 - <w>
  //               + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
  //           + ( 1 - mixt_frac ) * precip_frac_2
  //             * ( mu_w_2 - <w>
  //                 + corr_w_hm_2_n * sigma_w_2 * sigma_hm_2_n ) * mu_hm_2.
  // The overall covariance is provided, so the component correlation is solved
  // by setting corr_w_hm_1_n = corr_w_hm_2_n ( = corr_w_hm_n ).  The equation
  // is:
  // corr_w_hm_n
  // = ( <w'hm'>
  //     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
  //     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
  //   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1
  //       + ( 1 - mixt_frac ) * precip_frac_2
  //         * sigma_w_2 * sigma_hm_2_n * mu_hm_2 );
  // again, where corr_w_hm_1_n = corr_w_hm_2_n = corr_w_hm_n.  When either w
  // or hm is constant in one PDF component, but both w and hm vary in the
  // other PDF component, the equation for <w'hm'> is written as:
  // <w'hm'> = mixt_frac * precip_frac_1
  //           * ( mu_w_1 - <w>
  //               + corr_w_hm_1_n * sigma_w_1 * sigma_hm_1_n ) * mu_hm_1
  //           + ( 1 - mixt_frac ) * precip_frac_2
  //             * ( mu_w_2 - <w> ) * mu_hm_2.
  // In the above equation, either w or hm (or both) is (are) constant in PDF
  // component 2, but both w and hm vary in PDF component 1.  When both w and
  // hm vary in PDF component 2, but at least one of w or hm is constant in PDF
  // component 1, the equation is similar.  The above equation can be rewritten
  // to solve for corr_w_hm_1_n, such that:
  // corr_w_hm_1_n
  // = ( <w'hm'>
  //     - mixt_frac * precip_frac_1 * ( mu_w_1 - <w> ) * mu_hm_1
  //     - ( 1 - mixt_frac ) * precip_frac_2 * ( mu_w_2 - <w> ) * mu_hm_2 )
  //   / ( mixt_frac * precip_frac_1 * sigma_w_1 * sigma_hm_1_n * mu_hm_1 ).
  // Since either w or hm is constant in PDF component 2, corr_w_hm_2_n is
  // undefined.  When both w and hm vary in PDF component 2, but at least one
  // of w or hm is constant in PDF component 1, the equation is similar, but
  // is in terms of corr_w_hm_2_n, while corr_w_hm_1_n is undefined.  When
  // either w or hm is constant in both PDF components, the equation for
  // <w'hm'> is:
  // <w'hm'> = mixt_frac * precip_frac_1
  //           * ( mu_w_1 - <w> ) * mu_hm_1
  //           + ( 1 - mixt_frac ) * precip_frac_2
  //             * ( mu_w_2 - <w> ) * mu_hm_2.
  // When this is the case, both corr_w_hm_1_n and corr_w_hm_2_n are undefined.
  // References:
  //-----------------------------------------------------------------------

  // wm: Mean vertical velocity (overall), <w>             [m/s]
  // wphydrometp: Covariance of w and hm (overall), <w'hm'>  [m/s(hm un)]
  // mu_w_1: Mean of w (1st PDF component)                     [m/s]
  // mu_w_2: Mean of w (2nd PDF component)                     [m/s]
  // mu_hm_1: Mean of hm (1st PDF component) in-precip (ip)   [hm un]
  // mu_hm_2: Mean of hm (2nd PDF component) ip               [hm un]
  // sigma_w_1: Standard deviation of w (1st PDF component)       [m/s]
  // sigma_w_2: Standard deviation of w (2nd PDF component)       [m/s]
  // sigma_hm_1: Standard deviation of hm (1st PDF component) ip [hm un]
  // sigma_hm_2: Standard deviation of hm (2nd PDF component) ip [hm un]
  // sigma_hm_1_n: Standard deviation of ln hm (1st PDF component) ip  [-]
  // sigma_hm_2_n: Standard deviation of ln hm (2nd PDF component) ip  [-]
  // mixt_frac: Mixture fraction                                    [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)          [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)          [-]
  // hm_tol: Hydrometeor tolerance value                     [hm un]
  // corr_w_hm_1_n: Correlation of w and ln hm (1st PDF component) ip   [-]
  // corr_w_hm_2_n: Correlation of w and ln hm (2nd PDF component) ip   [-]

  real_t corr_w_hm_n{}; // Correlation of w and ln hm (both PDF components) ip   [-]

  // Calculate the PDF component correlation of vertical velocity, w, and the
  // natural logarithm of a hydrometeor, ln hm, in precipitation.
  if (((sigma_w_1 > w_tol) && (sigma_hm_1 > hm_tol) && (sigma_w_2 > w_tol) &&
       (sigma_hm_2 > hm_tol))) {
    // Both w and hm vary in both PDF components.
    // Calculate corr_w_hm_n (where corr_w_hm_1_n = corr_w_hm_2_n
    // = corr_w_hm_n).
    corr_w_hm_n =
        (((wphydrometp - (((mixt_frac * precip_frac_1) * (mu_w_1 - wm)) * mu_hm_1)) -
          ((((one - mixt_frac) * precip_frac_2) * (mu_w_2 - wm)) * mu_hm_2)) /
         (((((mixt_frac * precip_frac_1) * sigma_w_1) * sigma_hm_1_n) * mu_hm_1) +
          (((((one - mixt_frac) * precip_frac_2) * sigma_w_2) * sigma_hm_2_n) * mu_hm_2)));
    // Check that the PDF component correlations have reasonable values.
    if ((corr_w_hm_n > max_mag_correlation)) {
      corr_w_hm_n = max_mag_correlation;
    } else if ((corr_w_hm_n < -(max_mag_correlation))) {
      corr_w_hm_n = -(max_mag_correlation);
    }
    // The PDF component correlations between w and ln hm (in-precip) are
    // equal.
    corr_w_hm_1_n = corr_w_hm_n;
    corr_w_hm_2_n = corr_w_hm_n;
  } else if (((sigma_w_1 > w_tol) && (sigma_hm_1 > hm_tol))) {
    // Both w and hm vary in PDF component 1, but at least one of w and hm is
    // constant in PDF component 2.
    // Calculate the PDF component 1 correlation of w and ln hm (in-precip).
    corr_w_hm_1_n = (((wphydrometp - (((mixt_frac * precip_frac_1) * (mu_w_1 - wm)) * mu_hm_1)) -
                      ((((one - mixt_frac) * precip_frac_2) * (mu_w_2 - wm)) * mu_hm_2)) /
                     ((((mixt_frac * precip_frac_1) * sigma_w_1) * sigma_hm_1_n) * mu_hm_1));
    // Check that the PDF component 1 correlation has a reasonable value.
    if ((corr_w_hm_1_n > max_mag_correlation)) {
      corr_w_hm_1_n = max_mag_correlation;
    } else if ((corr_w_hm_1_n < -(max_mag_correlation))) {
      corr_w_hm_1_n = -(max_mag_correlation);
    }
    // The PDF component 2 correlation is undefined.
    corr_w_hm_2_n = zero;
  } else if (((sigma_w_2 > w_tol) && (sigma_hm_2 > hm_tol))) {
    // Both w and hm vary in PDF component 2, but at least one of w and hm is
    // constant in PDF component 1.
    // Calculate the PDF component 2 correlation of w and ln hm (in-precip).
    corr_w_hm_2_n =
        (((wphydrometp - (((mixt_frac * precip_frac_1) * (mu_w_1 - wm)) * mu_hm_1)) -
          ((((one - mixt_frac) * precip_frac_2) * (mu_w_2 - wm)) * mu_hm_2)) /
         (((((one - mixt_frac) * precip_frac_2) * sigma_w_2) * sigma_hm_2_n) * mu_hm_2));
    // Check that the PDF component 2 correlation has a reasonable value.
    if ((corr_w_hm_2_n > max_mag_correlation)) {
      corr_w_hm_2_n = max_mag_correlation;
    } else if ((corr_w_hm_2_n < -(max_mag_correlation))) {
      corr_w_hm_2_n = -(max_mag_correlation);
    }
    // The PDF component 1 correlation is undefined.
    corr_w_hm_1_n = zero;
  } else {
    // sigma_w_1 * sigma_hm_1 = 0 .and. sigma_w_2 * sigma_hm_2 = 0.
    // At least one of w and hm is constant in both PDF components.
    // The PDF component 1 and component 2 correlations are both undefined.
    corr_w_hm_1_n = zero;
    corr_w_hm_2_n = zero;
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void pdf_param_hm_stats(int nzt, int ngrdcol, int pdf_dim, int hydromet_dim, // In
                        const hm_metadata_type& hm_metadata,                 // In
                        Real3D hm_1, Real3D hm_2,                            // In
                        Real3D mu_x_1, Real3D mu_x_2,                        // In
                        Real3D sigma_x_1, Real3D sigma_x_2,                  // In
                        Real4D corr_array_1, Real4D corr_array_2,            // In
                        stats_type& stats                                    // InOut
) {
  // Description:
  // Record statistics for standard PDF parameters involving hydrometeors.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical levels
  // ngrdcol: Number of grid columns
  // pdf_dim: Number of variables in the correlation array
  // hydromet_dim: Number of hydrometeor species
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // hm_1: Mean of a precip. hydrometeor (1st PDF component)   [units vary]
  // hm_2: Mean of a precip. hydrometeor (2nd PDF component)   [units vary]
  // mu_x_1: Mean array of PDF vars. (1st PDF component)    [units vary]
  // mu_x_2: Mean array of PDF vars. (2nd PDF component)    [units vary]
  // sigma_x_1: Standard deviation array of PDF vars (comp. 1) [units vary]
  // sigma_x_2: Standard deviation array of PDF vars (comp. 2) [units vary]
  // corr_array_1: Correlation array of PDF vars. (comp. 1)             [-]
  // corr_array_2: Correlation array of PDF vars. (comp. 2)             [-]
  // stats: Native statistics collector.

  int ivar{};
  int jvar{};
  int hm_idx{};
  int hm_idx_ivar{};
  int hm_idx_jvar{};
  std::string hm_type{};
  std::string hmx_type{};
  std::string hmy_type{};
  std::string var_name{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_w{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi          = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta          = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w            = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn          = hm_metadata.iiPDF_Ncn;
  const auto hm_metadata_hydromet_list      = hm_metadata.hydromet_list;
  const auto hm_metadata_hydromet_list_host = to_host_space(hm_metadata_hydromet_list);

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // !! Output the statistics for hydrometeor PDF parameters.
  // Statistics
  if (stats.l_sample) {
    for (ivar = 0; ivar < hydromet_dim; ++ivar) {
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * ivar)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Mean of the precipitating hydrometeor in PDF component 1.
      var_name = setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name, Kokkos::subview(hm_1, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                            // InOut
      );
      // Mean of the precipitating hydrometeor in PDF component 2.
      var_name = setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name, Kokkos::subview(hm_2, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                            // InOut
      );
    }
    // Native loop: ivar = 1, hydromet_dim
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Mean of the precipitating hydrometeor (in-precip) in PDF component 1.
      var_name = "mu_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name, Kokkos::subview(mu_x_1, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                              // InOut
      );
      // Mean of the precipitating hydrometeor (in-precip) in PDF component 2.
      var_name = "mu_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name, Kokkos::subview(mu_x_2, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                              // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    // Mean of cloud nuclei concentration in PDF component 1.
    stats_update("mu_Ncn_1",                                                   // In
                 Kokkos::subview(mu_x_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                         // InOut
    );
    // Mean of cloud nuclei concentration in PDF component 2.
    stats_update("mu_Ncn_2",                                                   // In
                 Kokkos::subview(mu_x_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                         // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Standard deviation of the precipitating hydrometeor (in-precip) in PDF component 1.
      var_name = "sigma_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name, Kokkos::subview(sigma_x_1, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                                 // InOut
      );
      // Standard deviation of the precipitating hydrometeor (in-precip) in PDF component 2.
      var_name = "sigma_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name, Kokkos::subview(sigma_x_2, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                                 // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    // Standard deviation of cloud nuclei concentration in PDF component 1.
    stats_update("sigma_Ncn_1",                                                   // In
                 Kokkos::subview(sigma_x_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                            // InOut
    );
    // Standard deviation of cloud nuclei concentration in PDF component 2.
    stats_update("sigma_Ncn_2",                                                   // In
                 Kokkos::subview(sigma_x_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                            // InOut
    );
    // Correlations of w and chi/eta found in the correlation arrays.
    stats_update("corr_w_chi_1_ca",                                                           // In
                 Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_chi), // In
                 stats // InOut
    );
    stats_update("corr_w_chi_2_ca",                                                           // In
                 Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_chi), // In
                 stats // InOut
    );
    stats_update("corr_w_eta_1_ca",                                                           // In
                 Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_eta), // In
                 stats // InOut
    );
    stats_update("corr_w_eta_2_ca",                                                           // In
                 Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_w, iiPDF_eta), // In
                 stats // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of w and the precipitating hydrometeor in PDF component 1.
      var_name = "corr_w_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name,                                                               // In
                   Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_w), // In
                   stats                                                                   // InOut
      );
      // Correlation (in-precip) of w and the precipitating hydrometeor in PDF component 2.
      var_name = "corr_w_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name,                                                               // In
                   Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_w), // In
                   stats                                                                   // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    // Correlation of w and N_cn in PDF component 1.
    stats_update("corr_w_Ncn_1",                                                              // In
                 Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
                 stats // InOut
    );
    // Correlation of w and N_cn in PDF component 2.
    stats_update("corr_w_Ncn_2",                                                              // In
                 Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
                 stats // InOut
    );
    // Correlation of chi and eta in each PDF component found in the correlation arrays.
    stats_update(
        "corr_chi_eta_1_ca",                                                           // In
        Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_eta, iiPDF_chi), // In
        stats                                                                          // InOut
    );
    stats_update(
        "corr_chi_eta_2_ca",                                                           // In
        Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_eta, iiPDF_chi), // In
        stats                                                                          // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of chi (old s) and the precipitating
      // hydrometeor in PDF component 1.
      var_name = "corr_chi_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_chi), // In
                   stats // InOut
      );
      // Correlation (in-precip) of chi (old s) and the precipitating
      // hydrometeor in PDF component 2.
      var_name = "corr_chi_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_chi), // In
                   stats // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    // Correlation of chi (old s) and N_cn in PDF component 1.
    stats_update(
        "corr_chi_Ncn_1",                                                              // In
        Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi), // In
        stats                                                                          // InOut
    );
    // Correlation of chi (old s) and N_cn in PDF component 2.
    stats_update(
        "corr_chi_Ncn_2",                                                              // In
        Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_chi), // In
        stats                                                                          // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of eta (old t) and the precipitating
      // hydrometeor in PDF component 1.
      var_name = "corr_eta_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_eta), // In
                   stats // InOut
      );
      // Correlation (in-precip) of eta (old t) and the precipitating
      // hydrometeor in PDF component 2.
      var_name = "corr_eta_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_eta), // In
                   stats // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    // Correlation of eta (old t) and N_cn in PDF component 1.
    stats_update(
        "corr_eta_Ncn_1",                                                              // In
        Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta), // In
        stats                                                                          // InOut
    );
    // Correlation of eta (old t) and N_cn in PDF component 2.
    stats_update(
        "corr_eta_Ncn_2",                                                              // In
        Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_eta), // In
        stats                                                                          // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of N_cn and the precipitating hydrometeor in PDF component 1.
      var_name = "corr_Ncn_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_Ncn), // In
                   stats // InOut
      );
      // Correlation (in-precip) of N_cn and the precipitating hydrometeor in PDF component 2.
      var_name = "corr_Ncn_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_Ncn), // In
                   stats // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx_ivar = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hmx_type    = std::string(hm_metadata_hydromet_list_host.data() +
                                    ((hm_metadata_hydromet_list_host.extent(0) * hm_idx_ivar)),
                                hm_metadata_hydromet_list_host.extent(0));
      for (jvar = ((ivar + 1)); jvar < pdf_dim; ++jvar) {
        hm_idx_jvar = (pdf2hydromet_idx(jvar, hm_metadata) + 1) - 1;
        hmy_type    = std::string(hm_metadata_hydromet_list_host.data() +
                                      ((hm_metadata_hydromet_list_host.extent(0) * hm_idx_jvar)),
                                  hm_metadata_hydromet_list_host.extent(0));
        // Correlation (in-precip) of two different hydrometeors in PDF component 1.
        var_name = "corr_" + setup_clubb_pdf_params_detail::trim(hmx_type.substr(0, 2)) + "_" +
                   setup_clubb_pdf_params_detail::trim(hmy_type.substr(0, 2)) + "_1";
        stats_update(var_name,                                                            // In
                     Kokkos::subview(corr_array_1, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                     stats                                                                // InOut
        );
        // Correlation (in-precip) of two different hydrometeors in PDF component 2.
        var_name = "corr_" + setup_clubb_pdf_params_detail::trim(hmx_type.substr(0, 2)) + "_" +
                   setup_clubb_pdf_params_detail::trim(hmy_type.substr(0, 2)) + "_2";
        stats_update(var_name,                                                            // In
                     Kokkos::subview(corr_array_2, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                     stats                                                                // InOut
        );
      }
    }
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void pdf_param_ln_hm_stats(ScratchWorkspace& workspace, int nzt, int ngrdcol, int pdf_dim,
                           const hm_metadata_type& hm_metadata,    // In
                           Real3D mu_x_1_n, Real3D mu_x_2_n,       // In
                           Real3D sigma_x_1_n, Real3D sigma_x_2_n, // In
                           Real4D corr_array_1_n,                  // In
                           Real4D corr_array_2_n,                  // In
                           stats_type& stats                       // InOut
) {
  // Description:
  // Record statistics for normal space PDF parameters involving hydrometeors.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical levels
  // ngrdcol: Number of grid columns
  // pdf_dim: Number of variables in the correlation array
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // mu_x_1_n: Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // mu_x_2_n: Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // stats: Native statistics collector.

  auto workspace_scope = workspace.scope();
  Real2D mu_hm_1_n =
      workspace_scope.view<Real2D>("pdf_param_ln_hm_stats::mu_hm_1_n", ngrdcol,
                                   nzt); // Mean of ln hm (1st PDF component)    [units vary]
  Real2D mu_hm_2_n =
      workspace_scope.view<Real2D>("pdf_param_ln_hm_stats::mu_hm_2_n", ngrdcol,
                                   nzt); // Mean of ln hm (2nd PDF component)    [units vary]
  Real2D mu_Ncn_1_n =
      workspace_scope.view<Real2D>("pdf_param_ln_hm_stats::mu_Ncn_1_n", ngrdcol,
                                   nzt); // Mean of ln Ncn (1st PDF component)   [ln(num/kg)]
  Real2D mu_Ncn_2_n =
      workspace_scope.view<Real2D>("pdf_param_ln_hm_stats::mu_Ncn_2_n", ngrdcol,
                                   nzt); // Mean of ln Ncn (2nd PDF component)   [ln(num/kg)]
  int ivar{};
  int jvar{};
  int hm_idx{};
  int hm_idx_ivar{};
  int hm_idx_jvar{};
  std::string hm_type{};
  std::string hmx_type{};
  std::string hmy_type{};
  std::string var_name{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_w{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi          = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta          = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w            = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn          = hm_metadata.iiPDF_Ncn;
  const auto hm_metadata_hydromet_list      = hm_metadata.hydromet_list;
  const auto hm_metadata_hydromet_list_host = to_host_space(hm_metadata_hydromet_list);

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // !! Output the statistics for normal space hydrometeor PDF parameters.
  // Statistics
  if (stats.l_sample) {
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx   = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type  = std::string(hm_metadata_hydromet_list_host.data() +
                                 ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                             hm_metadata_hydromet_list_host.extent(0));
      var_name = "mu_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      if (var_on_stats_list(stats, var_name)) {
        {
          const int extent_0 = nzt;
          const int extent_1 = ngrdcol;
          Kokkos::parallel_for(
              "pdf_param_ln_hm_stats::cells_1", ((extent_0 * extent_1)),
              KOKKOS_LAMBDA(const int cell) {
                const int i = ((cell % extent_1));
                const int k = ((cell / extent_1));
                if ((mu_x_1_n(i, k, ivar) > -real_t(Kokkos::finite_max_v<float>))) {
                  mu_hm_1_n(i, k) = mu_x_1_n(i, k, ivar);
                } else {
                  // When hm_1 is 0 (or below tolerance value), mu_hm_1_n is -inf.
                  // Store -huge so output remains finite across precisions.
                  mu_hm_1_n(i, k) = -real_t(Kokkos::finite_max_v<float>);
                }
              });
        }
        stats_update(var_name, mu_hm_1_n, // In
                     stats                // InOut
        );
      }
      var_name = "mu_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      if (var_on_stats_list(stats, var_name)) {
        {
          const int extent_0 = nzt;
          const int extent_1 = ngrdcol;
          Kokkos::parallel_for(
              "pdf_param_ln_hm_stats::cells_2", ((extent_0 * extent_1)),
              KOKKOS_LAMBDA(const int cell) {
                const int i = ((cell % extent_1));
                const int k = ((cell / extent_1));
                if ((mu_x_2_n(i, k, ivar) > -real_t(Kokkos::finite_max_v<float>))) {
                  mu_hm_2_n(i, k) = mu_x_2_n(i, k, ivar);
                } else {
                  // When hm_2 is 0 (or below tolerance value), mu_hm_2_n is -inf.
                  // Store -huge so output remains finite across precisions.
                  mu_hm_2_n(i, k) = -real_t(Kokkos::finite_max_v<float>);
                }
              });
        }
        stats_update(var_name, mu_hm_2_n, // In
                     stats                // InOut
        );
      }
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    if (var_on_stats_list(stats, "mu_Ncn_1_n")) {
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "pdf_param_ln_hm_stats::cells_3", ((extent_0 * extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              if ((mu_x_1_n(i, k, iiPDF_Ncn) > -real_t(Kokkos::finite_max_v<float>))) {
                mu_Ncn_1_n(i, k) = mu_x_1_n(i, k, iiPDF_Ncn);
              } else {
                // When Ncnm is 0 (or below tolerance value), mu_Ncn_1_n is -inf.
                // Store -huge so output remains finite across precisions.
                mu_Ncn_1_n(i, k) = -real_t(Kokkos::finite_max_v<float>);
              }
            });
      }
      stats_update("mu_Ncn_1_n", mu_Ncn_1_n, // In
                   stats                     // InOut
      );
    }
    if (var_on_stats_list(stats, "mu_Ncn_2_n")) {
      {
        const int extent_0 = nzt;
        const int extent_1 = ngrdcol;
        Kokkos::parallel_for(
            "pdf_param_ln_hm_stats::cells_4", ((extent_0 * extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int i = ((cell % extent_1));
              const int k = ((cell / extent_1));
              if ((mu_x_2_n(i, k, iiPDF_Ncn) > -real_t(Kokkos::finite_max_v<float>))) {
                mu_Ncn_2_n(i, k) = mu_x_2_n(i, k, iiPDF_Ncn);
              } else {
                // When Ncnm is 0 (or below tolerance value), mu_Ncn_2_n is -inf.
                // Store -huge so output remains finite across precisions.
                mu_Ncn_2_n(i, k) = -real_t(Kokkos::finite_max_v<float>);
              }
            });
      }
      stats_update("mu_Ncn_2_n", mu_Ncn_2_n, // In
                   stats                     // InOut
      );
    }
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx   = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type  = std::string(hm_metadata_hydromet_list_host.data() +
                                 ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                             hm_metadata_hydromet_list_host.extent(0));
      var_name = "sigma_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      stats_update(var_name,                                                     // In
                   Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                         // InOut
      );
      var_name = "sigma_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      stats_update(var_name,                                                     // In
                   Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, ivar), // In
                   stats                                                         // InOut
      );
    }
    // Native loop: ivar = iiPDF_Ncn+1, pdf_dim, 1
    stats_update("sigma_Ncn_1_n",                                                   // In
                 Kokkos::subview(sigma_x_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                              // InOut
    );
    stats_update("sigma_Ncn_2_n",                                                   // In
                 Kokkos::subview(sigma_x_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn), // In
                 stats                                                              // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of w and ln hm in PDF component 1.
      var_name = "corr_w_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_w), // In
                   stats // InOut
      );
      // Correlation (in-precip) of w and ln hm in PDF component 2.
      var_name = "corr_w_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      stats_update(var_name,                                                                 // In
                   Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_w), // In
                   stats // InOut
      );
    }
    // Correlation of w and ln N_cn in PDF component 1.
    stats_update(
        "corr_w_Ncn_1_n",                                                              // In
        Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
        stats                                                                          // InOut
    );
    // Correlation of w and ln N_cn in PDF component 2.
    stats_update(
        "corr_w_Ncn_2_n",                                                              // In
        Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn, iiPDF_w), // In
        stats                                                                          // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of chi (old s) and ln hm in PDF component 1.
      var_name = "corr_chi_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_chi), // In
                   stats // InOut
      );
      // Correlation (in-precip) of chi (old s) and ln hm in PDF component 2.
      var_name = "corr_chi_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_chi), // In
                   stats // InOut
      );
    }
    // Correlation of chi (old s) and ln N_cn in PDF component 1.
    stats_update("corr_chi_Ncn_1_n", // In
                 Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                                 iiPDF_chi), // In
                 stats                       // InOut
    );
    // Correlation of chi (old s) and ln N_cn in PDF component 2.
    stats_update("corr_chi_Ncn_2_n", // In
                 Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                                 iiPDF_chi), // In
                 stats                       // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of eta (old t) and ln hm in PDF component 1.
      var_name = "corr_eta_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_eta), // In
                   stats // InOut
      );
      // Correlation (in-precip) of eta (old t) and ln hm in PDF component 2.
      var_name = "corr_eta_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_eta), // In
                   stats // InOut
      );
    }
    // Correlation of eta (old t) and ln N_cn in PDF component 1.
    stats_update("corr_eta_Ncn_1_n", // In
                 Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                                 iiPDF_eta), // In
                 stats                       // InOut
    );
    // Correlation of eta (old t) and ln N_cn in PDF component 2.
    stats_update("corr_eta_Ncn_2_n", // In
                 Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, iiPDF_Ncn,
                                 iiPDF_eta), // In
                 stats                       // InOut
    );
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx  = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hm_type = std::string(hm_metadata_hydromet_list_host.data() +
                                ((hm_metadata_hydromet_list_host.extent(0) * hm_idx)),
                            hm_metadata_hydromet_list_host.extent(0));
      // Correlation (in-precip) of ln N_cn and ln hm in PDF component 1.
      var_name = "corr_Ncn_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_1_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_Ncn), // In
                   stats // InOut
      );
      // Correlation (in-precip) of ln N_cn and ln hm in PDF component 2.
      var_name = "corr_Ncn_" + setup_clubb_pdf_params_detail::trim(hm_type.substr(0, 2)) + "_2_n";
      stats_update(var_name,                                                                   // In
                   Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, ivar, iiPDF_Ncn), // In
                   stats // InOut
      );
    }
    for (ivar = ((iiPDF_Ncn + 1)); ivar < pdf_dim; ++ivar) {
      hm_idx_ivar = (pdf2hydromet_idx(ivar, hm_metadata) + 1) - 1;
      hmx_type    = std::string(hm_metadata_hydromet_list_host.data() +
                                    ((hm_metadata_hydromet_list_host.extent(0) * hm_idx_ivar)),
                                hm_metadata_hydromet_list_host.extent(0));
      for (jvar = ((ivar + 1)); jvar < pdf_dim; ++jvar) {
        hm_idx_jvar = (pdf2hydromet_idx(jvar, hm_metadata) + 1) - 1;
        hmy_type    = std::string(hm_metadata_hydromet_list_host.data() +
                                      ((hm_metadata_hydromet_list_host.extent(0) * hm_idx_jvar)),
                                  hm_metadata_hydromet_list_host.extent(0));
        // Correlation (in-precip) of ln hmx and ln hmy in PDF component 1.
        var_name = "corr_" + setup_clubb_pdf_params_detail::trim(hmx_type.substr(0, 2)) + "_" +
                   setup_clubb_pdf_params_detail::trim(hmy_type.substr(0, 2)) + "_1_n";
        stats_update(var_name,                                                              // In
                     Kokkos::subview(corr_array_1_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                     stats                                                                  // InOut
        );
        // Correlation (in-precip) of ln hmx and ln hmy in PDF component 2.
        var_name = "corr_" + setup_clubb_pdf_params_detail::trim(hmx_type.substr(0, 2)) + "_" +
                   setup_clubb_pdf_params_detail::trim(hmy_type.substr(0, 2)) + "_2_n";
        stats_update(var_name,                                                              // In
                     Kokkos::subview(corr_array_2_n, Kokkos::ALL, Kokkos::ALL, jvar, ivar), // In
                     stats                                                                  // InOut
        );
      }
    }
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void pack_hydromet_pdf_params(int nzt, int ngrdcol,                                  // In
                              int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                              Real3D hm_1, Real3D hm_2, int pdf_dim, Real3D mu_x_1,  // In
                              Real3D mu_x_2, Real3D sigma_x_1, Real3D sigma_x_2,     // In
                              Real4D corr_array_1, Real4D corr_array_2,              // In
                              hydromet_pdf_parameter_collection& hydromet_pdf_params // Out
) {
  // Description:
  // Pack the standard means and variances involving hydrometeors, as well as a
  // few other variables, into the structure hydromet_pdf_params.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // hm_metadata: Hydrometeor names, tolerances, variance ratios and zero-based PDF slots.
  // hm_1: Mean of a precip. hydrometeor (1st PDF component)  [units vary]
  // hm_2: Mean of a precip. hydrometeor (2nd PDF component)  [units vary]
  // pdf_dim: Number of variables in the mean/stdev arrays
  // mu_x_1: Mean array of PDF vars. (1st PDF component)    [units vary]
  // mu_x_2: Mean array of PDF vars. (2nd PDF component)    [units vary]
  // sigma_x_1: Standard deviation array of PDF vars (comp. 1) [units vary]
  // sigma_x_2: Standard deviation array of PDF vars (comp. 2) [units vary]
  // corr_array_1: Correlation array of PDF vars. (comp. 1)    [-]
  // corr_array_2: Correlation array of PDF vars. (comp. 2)    [-]
  // hydromet_pdf_params: Hydrometeor PDF parameters        [units vary]

  int ivar{};
  int jvar{};
  int pdf_idx{};
  int iiPDF_chi{};
  int iiPDF_eta{};
  int iiPDF_w{};
  int iiPDF_Ncn{};
  // Capture component Views/scalars with zero-based PDF slots.
  // Owning UDTs remain on the host.
  const auto hm_metadata_iiPDF_chi              = hm_metadata.iiPDF_chi;
  const auto hm_metadata_iiPDF_eta              = hm_metadata.iiPDF_eta;
  const auto hm_metadata_iiPDF_w                = hm_metadata.iiPDF_w;
  const auto hm_metadata_iiPDF_Ncn              = hm_metadata.iiPDF_Ncn;
  const auto hydromet_pdf_params_hm_1           = hydromet_pdf_params.hm_1;
  const auto hydromet_pdf_params_hm_2           = hydromet_pdf_params.hm_2;
  const auto hydromet_pdf_params_mu_hm_1        = hydromet_pdf_params.mu_hm_1;
  const auto hydromet_pdf_params_mu_hm_2        = hydromet_pdf_params.mu_hm_2;
  const auto hydromet_pdf_params_sigma_hm_1     = hydromet_pdf_params.sigma_hm_1;
  const auto hydromet_pdf_params_sigma_hm_2     = hydromet_pdf_params.sigma_hm_2;
  const auto hydromet_pdf_params_corr_w_hm_1    = hydromet_pdf_params.corr_w_hm_1;
  const auto hydromet_pdf_params_corr_w_hm_2    = hydromet_pdf_params.corr_w_hm_2;
  const auto hydromet_pdf_params_corr_chi_hm_1  = hydromet_pdf_params.corr_chi_hm_1;
  const auto hydromet_pdf_params_corr_chi_hm_2  = hydromet_pdf_params.corr_chi_hm_2;
  const auto hydromet_pdf_params_corr_eta_hm_1  = hydromet_pdf_params.corr_eta_hm_1;
  const auto hydromet_pdf_params_corr_eta_hm_2  = hydromet_pdf_params.corr_eta_hm_2;
  const auto hydromet_pdf_params_corr_hmx_hmy_1 = hydromet_pdf_params.corr_hmx_hmy_1;
  const auto hydromet_pdf_params_corr_hmx_hmy_2 = hydromet_pdf_params.corr_hmx_hmy_2;
  const auto hydromet_pdf_params_mu_Ncn_1       = hydromet_pdf_params.mu_Ncn_1;
  const auto hydromet_pdf_params_mu_Ncn_2       = hydromet_pdf_params.mu_Ncn_2;
  const auto hydromet_pdf_params_sigma_Ncn_1    = hydromet_pdf_params.sigma_Ncn_1;
  const auto hydromet_pdf_params_sigma_Ncn_2    = hydromet_pdf_params.sigma_Ncn_2;

  iiPDF_chi = hm_metadata_iiPDF_chi;
  iiPDF_eta = hm_metadata_iiPDF_eta;
  iiPDF_w   = hm_metadata_iiPDF_w;
  iiPDF_Ncn = hm_metadata_iiPDF_Ncn;
  // Pack remaining means and standard deviations into hydromet_pdf_params.
  for (ivar = 0; ivar < hydromet_dim; ++ivar) {
    pdf_idx = (hydromet2pdf_idx(ivar, hm_metadata) + 1) - 1;
    // Mean of a hydrometeor (overall) in the 1st PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_1", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                = ((cell % slice_extent_0));
            const int slice_1                                = ((cell / slice_extent_0));
            hydromet_pdf_params_hm_1(slice_0, slice_1, ivar) = hm_1(slice_0, slice_1, ivar);
          });
    }
    // Mean of a hydrometeor (overall) in the 2nd PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_2", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                = ((cell % slice_extent_0));
            const int slice_1                                = ((cell / slice_extent_0));
            hydromet_pdf_params_hm_2(slice_0, slice_1, ivar) = hm_2(slice_0, slice_1, ivar);
          });
    }
    // Mean of a hydrometeor (in-precip) in the 1st PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_mu_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_mu_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_3", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                   = ((cell % slice_extent_0));
            const int slice_1                                   = ((cell / slice_extent_0));
            hydromet_pdf_params_mu_hm_1(slice_0, slice_1, ivar) = mu_x_1(slice_0, slice_1, pdf_idx);
          });
    }
    // Mean of a hydrometeor (in-precip) in the 2nd PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_mu_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_mu_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_4", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0                                   = ((cell % slice_extent_0));
            const int slice_1                                   = ((cell / slice_extent_0));
            hydromet_pdf_params_mu_hm_2(slice_0, slice_1, ivar) = mu_x_2(slice_0, slice_1, pdf_idx);
          });
    }
    // Standard deviation of a hydrometeor (in-precip) in the
    // 1st PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_sigma_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_sigma_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_5", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_sigma_hm_1(slice_0, slice_1, ivar) =
                sigma_x_1(slice_0, slice_1, pdf_idx);
          });
    }
    // Standard deviation of a hydrometeor (in-precip) in the
    // 2nd PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_sigma_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_sigma_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_6", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_sigma_hm_2(slice_0, slice_1, ivar) =
                sigma_x_2(slice_0, slice_1, pdf_idx);
          });
    }
    // Correlation (in-precip) of w and a hydrometeor in the 1st PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_w_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_w_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_7", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_w_hm_1(slice_0, slice_1, ivar) =
                corr_array_1(slice_0, slice_1, pdf_idx, iiPDF_w);
          });
    }
    // Correlation (in-precip) of w and a hydrometeor in the 2nd PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_w_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_w_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_8", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_w_hm_2(slice_0, slice_1, ivar) =
                corr_array_2(slice_0, slice_1, pdf_idx, iiPDF_w);
          });
    }
    // Correlation (in-precip) of chi and a hydrometeor in the 1st PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_chi_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_chi_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_9", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_chi_hm_1(slice_0, slice_1, ivar) =
                corr_array_1(slice_0, slice_1, pdf_idx, iiPDF_chi);
          });
    }
    // Correlation (in-precip) of chi and a hydrometeor in the 2nd PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_chi_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_chi_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_10", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_chi_hm_2(slice_0, slice_1, ivar) =
                corr_array_2(slice_0, slice_1, pdf_idx, iiPDF_chi);
          });
    }
    // Correlation (in-precip) of eta and a hydrometeor in the 1st PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_eta_hm_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_eta_hm_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_11", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_eta_hm_1(slice_0, slice_1, ivar) =
                corr_array_1(slice_0, slice_1, pdf_idx, iiPDF_eta);
          });
    }
    // Correlation (in-precip) of eta and a hydrometeor in the 2nd PDF
    // component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_eta_hm_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_eta_hm_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_12", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_eta_hm_2(slice_0, slice_1, ivar) =
                corr_array_2(slice_0, slice_1, pdf_idx, iiPDF_eta);
          });
    }
    // Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 1st
    // PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_13", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_hmx_hmy_1(slice_0, slice_1, ivar, ivar) = one;
          });
    }
    for (jvar = ((ivar + 1)); jvar < hydromet_dim; ++jvar) {
      {
        const int jpdf_idx = (hydromet2pdf_idx(jvar, hm_metadata) + 1) - 1;
        {
          const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(0);
          const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(1);
          Kokkos::parallel_for(
              "pack_hydromet_pdf_params::array_14", ((slice_extent_0 * slice_extent_1)),
              KOKKOS_LAMBDA(const int cell) {
                const int slice_0 = ((cell % slice_extent_0));
                const int slice_1 = ((cell / slice_extent_0));
                hydromet_pdf_params_corr_hmx_hmy_1(slice_0, slice_1, jvar, ivar) =
                    corr_array_1(slice_0, slice_1, jpdf_idx, pdf_idx);
              });
        }
      }
      {
        const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(0);
        const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_1.extent_int(1);
        Kokkos::parallel_for(
            "pack_hydromet_pdf_params::array_15", ((slice_extent_0 * slice_extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int slice_0 = ((cell % slice_extent_0));
              const int slice_1 = ((cell / slice_extent_0));
              hydromet_pdf_params_corr_hmx_hmy_1(slice_0, slice_1, ivar, jvar) =
                  hydromet_pdf_params_corr_hmx_hmy_1(slice_0, slice_1, jvar, ivar);
            });
      }
    }
    // Native loop: jvar = ivar+1, hydromet_dim, 1
    // Correlation (in-precip) of two hydrometeors, hmx and hmy, in the 2nd
    // PDF component.
    {
      const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(0);
      const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(1);
      Kokkos::parallel_for(
          "pack_hydromet_pdf_params::array_16", ((slice_extent_0 * slice_extent_1)),
          KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = ((cell % slice_extent_0));
            const int slice_1 = ((cell / slice_extent_0));
            hydromet_pdf_params_corr_hmx_hmy_2(slice_0, slice_1, ivar, ivar) = one;
          });
    }
    for (jvar = ((ivar + 1)); jvar < hydromet_dim; ++jvar) {
      {
        const int jpdf_idx = (hydromet2pdf_idx(jvar, hm_metadata) + 1) - 1;
        {
          const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(0);
          const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(1);
          Kokkos::parallel_for(
              "pack_hydromet_pdf_params::array_17", ((slice_extent_0 * slice_extent_1)),
              KOKKOS_LAMBDA(const int cell) {
                const int slice_0 = ((cell % slice_extent_0));
                const int slice_1 = ((cell / slice_extent_0));
                hydromet_pdf_params_corr_hmx_hmy_2(slice_0, slice_1, jvar, ivar) =
                    corr_array_2(slice_0, slice_1, jpdf_idx, pdf_idx);
              });
        }
      }
      {
        const int slice_extent_0 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(0);
        const int slice_extent_1 = hydromet_pdf_params_corr_hmx_hmy_2.extent_int(1);
        Kokkos::parallel_for(
            "pack_hydromet_pdf_params::array_18", ((slice_extent_0 * slice_extent_1)),
            KOKKOS_LAMBDA(const int cell) {
              const int slice_0 = ((cell % slice_extent_0));
              const int slice_1 = ((cell / slice_extent_0));
              hydromet_pdf_params_corr_hmx_hmy_2(slice_0, slice_1, ivar, jvar) =
                  hydromet_pdf_params_corr_hmx_hmy_2(slice_0, slice_1, jvar, ivar);
            });
      }
    }
    // Native loop: jvar = ivar+1, hydromet_dim, 1
  }
  // Native loop: ivar = 1, hydromet_dim, 1
  // Mean of Ncn (overall) in the 1st PDF component.
  {
    const int slice_extent_0 = hydromet_pdf_params_mu_Ncn_1.extent_int(0);
    const int slice_extent_1 = hydromet_pdf_params_mu_Ncn_1.extent_int(1);
    Kokkos::parallel_for(
        "pack_hydromet_pdf_params::array_19", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                              = ((cell % slice_extent_0));
          const int slice_1                              = ((cell / slice_extent_0));
          hydromet_pdf_params_mu_Ncn_1(slice_0, slice_1) = mu_x_1(slice_0, slice_1, iiPDF_Ncn);
        });
  }
  // Mean of Ncn (overall) in the 2nd PDF component.
  {
    const int slice_extent_0 = hydromet_pdf_params_mu_Ncn_2.extent_int(0);
    const int slice_extent_1 = hydromet_pdf_params_mu_Ncn_2.extent_int(1);
    Kokkos::parallel_for(
        "pack_hydromet_pdf_params::array_20", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0                              = ((cell % slice_extent_0));
          const int slice_1                              = ((cell / slice_extent_0));
          hydromet_pdf_params_mu_Ncn_2(slice_0, slice_1) = mu_x_2(slice_0, slice_1, iiPDF_Ncn);
        });
  }
  // Standard deviation of Ncn (overall) in the 1st PDF component.
  {
    const int slice_extent_0 = hydromet_pdf_params_sigma_Ncn_1.extent_int(0);
    const int slice_extent_1 = hydromet_pdf_params_sigma_Ncn_1.extent_int(1);
    Kokkos::parallel_for(
        "pack_hydromet_pdf_params::array_21", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          hydromet_pdf_params_sigma_Ncn_1(slice_0, slice_1) =
              sigma_x_1(slice_0, slice_1, iiPDF_Ncn);
        });
  }
  // Standard deviation of Ncn (overall) in the 2nd PDF component.
  {
    const int slice_extent_0 = hydromet_pdf_params_sigma_Ncn_2.extent_int(0);
    const int slice_extent_1 = hydromet_pdf_params_sigma_Ncn_2.extent_int(1);
    Kokkos::parallel_for(
        "pack_hydromet_pdf_params::array_22", ((slice_extent_0 * slice_extent_1)),
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = ((cell % slice_extent_0));
          const int slice_1 = ((cell / slice_extent_0));
          hydromet_pdf_params_sigma_Ncn_2(slice_0, slice_1) =
              sigma_x_2(slice_0, slice_1, iiPDF_Ncn);
        });
  }
  return;
}
} // namespace
// =============================================================================
namespace {
KOKKOS_INLINE_FUNCTION real_t compute_rtp2_from_chi(real_t sigma_chi_1, real_t sigma_chi_2, // In
                                                    real_t sigma_eta_1, real_t sigma_eta_2, // In
                                                    real_t rt_1, real_t rt_2,               // In
                                                    real_t crt_1, real_t crt_2,             // In
                                                    real_t mixt_frac,                       // In
                                                    real_t corr_chi_eta_1,
                                                    real_t corr_chi_eta_2 // In
) {
  // Description:
  // Compute the variance of rt from the distribution of chi and eta. The
  // resulting variance will be consistent with CLUBB's extended PDF
  // involving chi and eta, including if l_fix_w_chi_eta_correlations = .true..
  // References:
  //   None
  //-----------------------------------------------------------------------

  // sigma_chi_1: Standard deviation of chi (1st PDF comp.) [kg/kg]
  // sigma_chi_2: Standard deviation of chi (2nd PDF comp.) [kg/kg]
  // sigma_eta_1: Standard deviation of eta (1st PDF comp.) [kg/kg]
  // sigma_eta_2: Standard deviation of eta (2nd PDF comp.) [kg/kg]
  // rt_1: Mean of rt (1st PDF component)            [kg/kg]
  // rt_2: Mean of rt (2nd PDF component)            [kg/kg]
  // crt_1: Coef. of r_t in chi/eta eqns. (1st comp.) [-]
  // crt_2: Coef. of r_t in chi/eta eqns. (2nd comp.) [-]
  // mixt_frac: Weight of 1st gaussian PDF component      [-]
  // corr_chi_eta_1: Correlation of chi and eta in 1st PDF component [-]
  // corr_chi_eta_2: Correlation of chi and eta in 2nd PDF component [-]

  real_t rtp2_zt_from_chi{}; // Grid-box rt variance on thermo levels [(kg/kg)^2]. // Grid-box
                             // variance of rtp2 on thermo. levels  [kg/kg]
  real_t varnce_rt_1_zt_from_chi{};
  real_t varnce_rt_2_zt_from_chi{};
  real_t rtm{};                 // Mean of rt (overall) [kg/kg]
  real_t sigma_rt_1_from_chi{}; // Standard deviation of rt (1st PDF comp.)  [kg/kg]
  real_t sigma_rt_2_from_chi{}; // Standard deviation of rt (2nd PDF comp.)  [kg/kg]

  varnce_rt_1_zt_from_chi = (((((corr_chi_eta_1 * sigma_chi_1) * sigma_eta_1) +
                               (one_half * (sigma_chi_1 * sigma_chi_1))) +
                              (one_half * (sigma_eta_1 * sigma_eta_1))) /
                             (two * (crt_1 * crt_1)));
  varnce_rt_2_zt_from_chi = (((((corr_chi_eta_2 * sigma_chi_2) * sigma_eta_2) +
                               (one_half * (sigma_chi_2 * sigma_chi_2))) +
                              (one_half * (sigma_eta_2 * sigma_eta_2))) /
                             (two * (crt_2 * crt_2)));
  rtm                     = ((mixt_frac * rt_1) + ((one - mixt_frac) * rt_2));
  sigma_rt_1_from_chi     = Kokkos::sqrt(varnce_rt_1_zt_from_chi);
  sigma_rt_2_from_chi     = Kokkos::sqrt(varnce_rt_2_zt_from_chi);
  rtp2_zt_from_chi        = compute_variance_binormal(rtm, rt_1, rt_2, sigma_rt_1_from_chi,
                                                      sigma_rt_2_from_chi, mixt_frac);
  return rtp2_zt_from_chi;
}
} // namespace
} // namespace clubb_kokkos
