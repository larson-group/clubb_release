// Set extended CLUBB PDF moments, correlations and hydrometeor parameters.
//
// setup_pdf_parameters_api coordinates precipitation fractions, moments,
// normal-space transformations, matrix factors and Native-backed stats.
// Six public interfaces and twelve private helpers mirror the Native order.
// Fields use column-first LayoutLeft Views. Metadata indices stay one-based,
// while tunable parameter indices use their existing zero-based C++ owner.
// Scratch belongs to each call; optional hydrometeor outputs use the existing
// collection owner. Native diagnostic definedness and component call choices
// are documented at their source boundaries.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
#include "corr_varnce_module.hpp"
#include "pdf_parameter_module.hpp"
#include "hydromet_pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"
namespace clubb_kokkos {
// gr: Grid structure; physical bounds retain Native one-based indices.
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
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component PDF.
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
);
// nzt: Number of model thermodynamic vertical grid levels
// ngrdcol: Number of grid columns
// hydromet_dim: Number of hydrometeor species
// hydromet: Mean of a hydrometeor (overall)         [hm units]
// hydrometp2_zt: Variance of a hydrometeor (overall)     [(hm units)^2]
// hm_metadata: Hydrometeor names, tolerances, variance ratios and Native one-based PDF indices.
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
);
// nzt: Number of thermodynamic vertical levels
// pdf_dim: Number of variables in the corr/mean/stdev arrays
// ngrdcol: Number of grid columns
// hydromet_dim: Number of hydrometeor species
// hm_metadata: Hydrometeor names, tolerances, variance ratios and Native one-based PDF indices.
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
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component PDF.
// l_calc_w_corr: Calculate the correlations between w and the hydrometeors
// l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta)
// corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
// corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
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
);
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
);
// nzt: Number of model thermodynamic vertical grid levels
// ngrdcol: Number of grid columns
// hydromet_dim: Number of hydrometeor species
// hm_metadata: Hydrometeor names, tolerances, variance ratios and Native one-based PDF indices.
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
);
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// pdf_dim: Number of variables
// hm_metadata: Hydrometeor names, tolerances, variance ratios and Native one-based PDF indices.
// sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
// sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
// sigma2_on_mu2_ip_1: Ratio array sigma_hm_1^2/mu_hm_1^2             [-]
// sigma2_on_mu2_ip_2: Ratio array sigma_hm_2^2/mu_hm_2^2             [-]
// corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
// corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
// corr_array_1: Correlation array of PDF vars. (comp. 1)             [-]
// corr_array_2: Correlation array of PDF vars. (comp. 2)             [-]
void denorm_transform_corr(int nzt, int ngrdcol, int pdf_dim,
                           const hm_metadata_type& hm_metadata,                  // In
                           Real3D sigma_x_1_n, Real3D sigma_x_2_n,               // In
                           Real3D sigma2_on_mu2_ip_1, Real3D sigma2_on_mu2_ip_2, // In
                           Real4D corr_array_1_n,                                // In
                           Real4D corr_array_2_n,                                // In
                           Real4D corr_array_1, Real4D corr_array_2              // Out
);
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15, class Field16, class Field17>
inline void
setup_pdf_parameters_api(ScratchWorkspace& workspace,                                   // InOut
                         const grid& gr, int nzm, int nzt, int ngrdcol, int pdf_dim,    // In
                         int hydromet_dim,                                              // In
                         Field0 Nc_in_cloud, Field1 cloud_frac, Field2 Kh_zm,           // In
                         Field3 ice_supersat_frac, Field4 hydromet, Field5 wphydrometp, // In
                         Field6 corr_array_n_cloud, Field7 corr_array_n_below,          // In
                         const hm_metadata_type& hm_metadata,                           // In
                         const pdf_parameter& pdf_params,                               // In
                         Field8 clubb_params,                                           // In
                         int iiPDF_type,                                                // In
                         bool l_use_precip_frac,                                        // In
                         bool l_diagnose_correlations,                                  // In
                         bool l_calc_w_corr,                                            // In
                         bool l_const_Nc_in_cloud,                                      // In
                         bool l_fix_w_chi_eta_correlations,                             // In
                         err_info_type& err_info,                                       // InOut
                         Field9 hydrometp2,                                             // Out
                         Field10 mu_x_1_n, Field11 mu_x_2_n,                            // Out
                         Field12 sigma_x_1_n, Field13 sigma_x_2_n,                      // Out
                         Field14 corr_array_1_n, Field15 corr_array_2_n,                // Out
                         Field16 corr_cholesky_mtx_1, Field17 corr_cholesky_mtx_2,      // Out
                         precipitation_fractions& precip_fracs,                         // InOut
                         hydromet_pdf_parameter_collection* hydromet_pdf_params,        // Out
                         stats_type& stats                                              // InOut
) {
  auto workspace_scope                = workspace.scope();
  const auto Nc_in_cloud_exec         = to_execution_space(workspace, Nc_in_cloud);
  const auto cloud_frac_exec          = to_execution_space(workspace, cloud_frac);
  const auto Kh_zm_exec               = to_execution_space(workspace, Kh_zm);
  const auto ice_supersat_frac_exec   = to_execution_space(workspace, ice_supersat_frac);
  const auto hydromet_exec            = to_execution_space(workspace, hydromet);
  const auto wphydrometp_exec         = to_execution_space(workspace, wphydrometp);
  const auto corr_array_n_cloud_exec  = to_execution_space(workspace, corr_array_n_cloud);
  const auto corr_array_n_below_exec  = to_execution_space(workspace, corr_array_n_below);
  const auto clubb_params_exec        = to_execution_space(workspace, clubb_params);
  const auto hydrometp2_exec          = to_execution_space(workspace, hydrometp2, false);
  const auto mu_x_1_n_exec            = to_execution_space(workspace, mu_x_1_n, false);
  const auto mu_x_2_n_exec            = to_execution_space(workspace, mu_x_2_n, false);
  const auto sigma_x_1_n_exec         = to_execution_space(workspace, sigma_x_1_n, false);
  const auto sigma_x_2_n_exec         = to_execution_space(workspace, sigma_x_2_n, false);
  const auto corr_array_1_n_exec      = to_execution_space(workspace, corr_array_1_n, false);
  const auto corr_array_2_n_exec      = to_execution_space(workspace, corr_array_2_n, false);
  const auto corr_cholesky_mtx_1_exec = to_execution_space(workspace, corr_cholesky_mtx_1, false);
  const auto corr_cholesky_mtx_2_exec = to_execution_space(workspace, corr_cholesky_mtx_2, false);
  setup_pdf_parameters_api(
      workspace, gr, nzm, nzt, ngrdcol, pdf_dim, hydromet_dim, Nc_in_cloud_exec, cloud_frac_exec,
      Kh_zm_exec, ice_supersat_frac_exec, hydromet_exec, wphydrometp_exec, corr_array_n_cloud_exec,
      corr_array_n_below_exec, hm_metadata, pdf_params, clubb_params_exec, iiPDF_type,
      l_use_precip_frac, l_diagnose_correlations, l_calc_w_corr, l_const_Nc_in_cloud,
      l_fix_w_chi_eta_correlations, err_info, hydrometp2_exec, mu_x_1_n_exec, mu_x_2_n_exec,
      sigma_x_1_n_exec, sigma_x_2_n_exec, corr_array_1_n_exec, corr_array_2_n_exec,
      corr_cholesky_mtx_1_exec, corr_cholesky_mtx_2_exec, precip_fracs, hydromet_pdf_params, stats);
  from_execution_space(workspace, hydrometp2, hydrometp2_exec);
  from_execution_space(workspace, mu_x_1_n, mu_x_1_n_exec);
  from_execution_space(workspace, mu_x_2_n, mu_x_2_n_exec);
  from_execution_space(workspace, sigma_x_1_n, sigma_x_1_n_exec);
  from_execution_space(workspace, sigma_x_2_n, sigma_x_2_n_exec);
  from_execution_space(workspace, corr_array_1_n, corr_array_1_n_exec);
  from_execution_space(workspace, corr_array_2_n, corr_array_2_n_exec);
  from_execution_space(workspace, corr_cholesky_mtx_1, corr_cholesky_mtx_1_exec);
  from_execution_space(workspace, corr_cholesky_mtx_2, corr_cholesky_mtx_2_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15, class Field16, class Field17,
          class Field18, class Field19, class Field20, class Field21, class Field22, class Field23,
          class Field24, class Field25, class Field26, class Field27, class Field28, class Field29>
inline void compute_mean_stdev(int nzt, int ngrdcol, int hydromet_dim,   // In
                               Field0 hydromet, Field1 hydrometp2_zt,    // In
                               const hm_metadata_type& hm_metadata,      // In
                               Field2 Ncnm, Field3 mixt_frac,            // In
                               Field4 precip_frac,                       // In
                               Field5 precip_frac_1,                     // In
                               Field6 precip_frac_2,                     // In
                               Field7 precip_frac_tol,                   // In
                               Field8 w_1, Field9 w_2,                   // In
                               Field10 stdev_w_1, Field11 stdev_w_2,     // In
                               Field12 chi_1, Field13 chi_2,             // In
                               Field14 stdev_chi_1, Field15 stdev_chi_2, // In
                               Field16 stdev_eta_1, Field17 stdev_eta_2, // In
                               Field18 thl_1, Field19 thl_2,             // In
                               int pdf_dim,                              // In
                               Field20 omicron, Field21 zeta_vrnce_rat,  // In
                               bool l_const_Nc_in_cloud,                 // In
                               Field22 mu_x_1, Field23 mu_x_2,           // Out
                               Field24 sigma_x_1, Field25 sigma_x_2,     // Out
                               Field26 hm_1, Field27 hm_2,               // Out
                               Field28 sigma_hm_1_sqd_on_mu_hm_1_sqd,    // Out
                               Field29 sigma_hm_2_sqd_on_mu_hm_2_sqd     // Out
) {
  const auto hydromet_exec        = to_execution_space(hydromet);
  const auto hydrometp2_zt_exec   = to_execution_space(hydrometp2_zt);
  const auto Ncnm_exec            = to_execution_space(Ncnm);
  const auto mixt_frac_exec       = to_execution_space(mixt_frac);
  const auto precip_frac_exec     = to_execution_space(precip_frac);
  const auto precip_frac_1_exec   = to_execution_space(precip_frac_1);
  const auto precip_frac_2_exec   = to_execution_space(precip_frac_2);
  const auto precip_frac_tol_exec = to_execution_space(precip_frac_tol);
  const auto w_1_exec             = to_execution_space(w_1);
  const auto w_2_exec             = to_execution_space(w_2);
  const auto stdev_w_1_exec       = to_execution_space(stdev_w_1);
  const auto stdev_w_2_exec       = to_execution_space(stdev_w_2);
  const auto chi_1_exec           = to_execution_space(chi_1);
  const auto chi_2_exec           = to_execution_space(chi_2);
  const auto stdev_chi_1_exec     = to_execution_space(stdev_chi_1);
  const auto stdev_chi_2_exec     = to_execution_space(stdev_chi_2);
  const auto stdev_eta_1_exec     = to_execution_space(stdev_eta_1);
  const auto stdev_eta_2_exec     = to_execution_space(stdev_eta_2);
  const auto thl_1_exec           = to_execution_space(thl_1);
  const auto thl_2_exec           = to_execution_space(thl_2);
  const auto omicron_exec         = to_execution_space(omicron);
  const auto zeta_vrnce_rat_exec  = to_execution_space(zeta_vrnce_rat);
  const auto mu_x_1_exec          = to_execution_space(mu_x_1, false);
  const auto mu_x_2_exec          = to_execution_space(mu_x_2, false);
  const auto sigma_x_1_exec       = to_execution_space(sigma_x_1, false);
  const auto sigma_x_2_exec       = to_execution_space(sigma_x_2, false);
  const auto hm_1_exec            = to_execution_space(hm_1, false);
  const auto hm_2_exec            = to_execution_space(hm_2, false);
  const auto sigma_hm_1_sqd_on_mu_hm_1_sqd_exec =
      to_execution_space(sigma_hm_1_sqd_on_mu_hm_1_sqd, false);
  const auto sigma_hm_2_sqd_on_mu_hm_2_sqd_exec =
      to_execution_space(sigma_hm_2_sqd_on_mu_hm_2_sqd, false);
  compute_mean_stdev(nzt, ngrdcol, hydromet_dim, hydromet_exec, hydrometp2_zt_exec, hm_metadata,
                     Ncnm_exec, mixt_frac_exec, precip_frac_exec, precip_frac_1_exec,
                     precip_frac_2_exec, precip_frac_tol_exec, w_1_exec, w_2_exec, stdev_w_1_exec,
                     stdev_w_2_exec, chi_1_exec, chi_2_exec, stdev_chi_1_exec, stdev_chi_2_exec,
                     stdev_eta_1_exec, stdev_eta_2_exec, thl_1_exec, thl_2_exec, pdf_dim,
                     omicron_exec, zeta_vrnce_rat_exec, l_const_Nc_in_cloud, mu_x_1_exec,
                     mu_x_2_exec, sigma_x_1_exec, sigma_x_2_exec, hm_1_exec, hm_2_exec,
                     sigma_hm_1_sqd_on_mu_hm_1_sqd_exec, sigma_hm_2_sqd_on_mu_hm_2_sqd_exec);
  from_execution_space(mu_x_1, mu_x_1_exec);
  from_execution_space(mu_x_2, mu_x_2_exec);
  from_execution_space(sigma_x_1, sigma_x_1_exec);
  from_execution_space(sigma_x_2, sigma_x_2_exec);
  from_execution_space(hm_1, hm_1_exec);
  from_execution_space(hm_2, hm_2_exec);
  from_execution_space(sigma_hm_1_sqd_on_mu_hm_1_sqd, sigma_hm_1_sqd_on_mu_hm_1_sqd_exec);
  from_execution_space(sigma_hm_2_sqd_on_mu_hm_2_sqd, sigma_hm_2_sqd_on_mu_hm_2_sqd_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15, class Field16, class Field17>
inline void comp_corr_norm(ScratchWorkspace& workspace,                         // InOut
                           int nzt, int pdf_dim, int ngrdcol, int hydromet_dim, // In
                           const hm_metadata_type& hm_metadata,                 // In
                           Field0 wm_zt, Field1 rc_1, Field2 rc_2,              // In
                           Field3 mixt_frac,                                    // In
                           Field4 precip_frac_1,                                // In
                           Field5 precip_frac_2,                                // In
                           Field6 wpNcnp_zt, Field7 wphydrometp_zt,             // In
                           Field8 mu_x_1, Field9 mu_x_2, Field10 sigma_x_1, Field11 sigma_x_2, // In
                           Field12 sigma_x_1_n, Field13 sigma_x_2_n,                           // In
                           Field14 corr_array_n_cloud, Field15 corr_array_n_below,             // In
                           const pdf_parameter& pdf_params,                                    // In
                           int iiPDF_type,                                                     // In
                           bool l_calc_w_corr,                                                 // In
                           bool l_fix_w_chi_eta_correlations,                                  // In
                           Field16 corr_array_1_n, Field17 corr_array_2_n // Out
) {
  auto workspace_scope               = workspace.scope();
  const auto wm_zt_exec              = to_execution_space(workspace, wm_zt);
  const auto rc_1_exec               = to_execution_space(workspace, rc_1);
  const auto rc_2_exec               = to_execution_space(workspace, rc_2);
  const auto mixt_frac_exec          = to_execution_space(workspace, mixt_frac);
  const auto precip_frac_1_exec      = to_execution_space(workspace, precip_frac_1);
  const auto precip_frac_2_exec      = to_execution_space(workspace, precip_frac_2);
  const auto wpNcnp_zt_exec          = to_execution_space(workspace, wpNcnp_zt);
  const auto wphydrometp_zt_exec     = to_execution_space(workspace, wphydrometp_zt);
  const auto mu_x_1_exec             = to_execution_space(workspace, mu_x_1);
  const auto mu_x_2_exec             = to_execution_space(workspace, mu_x_2);
  const auto sigma_x_1_exec          = to_execution_space(workspace, sigma_x_1);
  const auto sigma_x_2_exec          = to_execution_space(workspace, sigma_x_2);
  const auto sigma_x_1_n_exec        = to_execution_space(workspace, sigma_x_1_n);
  const auto sigma_x_2_n_exec        = to_execution_space(workspace, sigma_x_2_n);
  const auto corr_array_n_cloud_exec = to_execution_space(workspace, corr_array_n_cloud);
  const auto corr_array_n_below_exec = to_execution_space(workspace, corr_array_n_below);
  const auto corr_array_1_n_exec     = to_execution_space(workspace, corr_array_1_n, false);
  const auto corr_array_2_n_exec     = to_execution_space(workspace, corr_array_2_n, false);
  comp_corr_norm(workspace, nzt, pdf_dim, ngrdcol, hydromet_dim, hm_metadata, wm_zt_exec, rc_1_exec,
                 rc_2_exec, mixt_frac_exec, precip_frac_1_exec, precip_frac_2_exec, wpNcnp_zt_exec,
                 wphydrometp_zt_exec, mu_x_1_exec, mu_x_2_exec, sigma_x_1_exec, sigma_x_2_exec,
                 sigma_x_1_n_exec, sigma_x_2_n_exec, corr_array_n_cloud_exec,
                 corr_array_n_below_exec, pdf_params, iiPDF_type, l_calc_w_corr,
                 l_fix_w_chi_eta_correlations, corr_array_1_n_exec, corr_array_2_n_exec);
  from_execution_space(workspace, corr_array_1_n, corr_array_1_n_exec);
  from_execution_space(workspace, corr_array_2_n, corr_array_2_n_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15, class Field16, class Field17,
          class Field18>
inline void calc_comp_mu_sigma_hm(int nzt, int ngrdcol,                       // In
                                  Field0 hmm, Field1 hmp2,                    // In
                                  real_t hmp2_ip_on_hmm2_ip,                  // In
                                  Field2 mixt_frac, Field3 precip_frac,       // In
                                  Field4 precip_frac_1, Field5 precip_frac_2, // In
                                  real_t hm_tol, Field6 precip_frac_tol,      // In
                                  Field7 mu_thl_1, Field8 mu_thl_2,           // In
                                  Field9 omicron, Field10 zeta_vrnce_rat_in,  // In
                                  Field11 mu_hm_1, Field12 mu_hm_2,           // Out
                                  Field13 sigma_hm_1, Field14 sigma_hm_2,     // Out
                                  Field15 hm_1, Field16 hm_2,                 // Out
                                  Field17 sigma_hm_1_sqd_on_mu_hm_1_sqd,      // Out
                                  Field18 sigma_hm_2_sqd_on_mu_hm_2_sqd       // Out
) {
  const auto hmm_exec               = to_execution_space(hmm);
  const auto hmp2_exec              = to_execution_space(hmp2);
  const auto mixt_frac_exec         = to_execution_space(mixt_frac);
  const auto precip_frac_exec       = to_execution_space(precip_frac);
  const auto precip_frac_1_exec     = to_execution_space(precip_frac_1);
  const auto precip_frac_2_exec     = to_execution_space(precip_frac_2);
  const auto precip_frac_tol_exec   = to_execution_space(precip_frac_tol);
  const auto mu_thl_1_exec          = to_execution_space(mu_thl_1);
  const auto mu_thl_2_exec          = to_execution_space(mu_thl_2);
  const auto omicron_exec           = to_execution_space(omicron);
  const auto zeta_vrnce_rat_in_exec = to_execution_space(zeta_vrnce_rat_in);
  const auto mu_hm_1_exec           = to_execution_space(mu_hm_1, false);
  const auto mu_hm_2_exec           = to_execution_space(mu_hm_2, false);
  const auto sigma_hm_1_exec        = to_execution_space(sigma_hm_1, false);
  const auto sigma_hm_2_exec        = to_execution_space(sigma_hm_2, false);
  const auto hm_1_exec              = to_execution_space(hm_1, false);
  const auto hm_2_exec              = to_execution_space(hm_2, false);
  const auto sigma_hm_1_sqd_on_mu_hm_1_sqd_exec =
      to_execution_space(sigma_hm_1_sqd_on_mu_hm_1_sqd, false);
  const auto sigma_hm_2_sqd_on_mu_hm_2_sqd_exec =
      to_execution_space(sigma_hm_2_sqd_on_mu_hm_2_sqd, false);
  calc_comp_mu_sigma_hm(nzt, ngrdcol, hmm_exec, hmp2_exec, hmp2_ip_on_hmm2_ip, mixt_frac_exec,
                        precip_frac_exec, precip_frac_1_exec, precip_frac_2_exec, hm_tol,
                        precip_frac_tol_exec, mu_thl_1_exec, mu_thl_2_exec, omicron_exec,
                        zeta_vrnce_rat_in_exec, mu_hm_1_exec, mu_hm_2_exec, sigma_hm_1_exec,
                        sigma_hm_2_exec, hm_1_exec, hm_2_exec, sigma_hm_1_sqd_on_mu_hm_1_sqd_exec,
                        sigma_hm_2_sqd_on_mu_hm_2_sqd_exec);
  from_execution_space(mu_hm_1, mu_hm_1_exec);
  from_execution_space(mu_hm_2, mu_hm_2_exec);
  from_execution_space(sigma_hm_1, sigma_hm_1_exec);
  from_execution_space(sigma_hm_2, sigma_hm_2_exec);
  from_execution_space(hm_1, hm_1_exec);
  from_execution_space(hm_2, hm_2_exec);
  from_execution_space(sigma_hm_1_sqd_on_mu_hm_1_sqd, sigma_hm_1_sqd_on_mu_hm_1_sqd_exec);
  from_execution_space(sigma_hm_2_sqd_on_mu_hm_2_sqd, sigma_hm_2_sqd_on_mu_hm_2_sqd_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12>
inline void norm_transform_mean_stdev(int nzt, int ngrdcol,                                  // In
                                      int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                                      Field0 hm_1, Field1 hm_2,                              // In
                                      Field2 Ncnm, int pdf_dim,                              // In
                                      Field3 mu_x_1, Field4 mu_x_2,                          // In
                                      Field5 sigma_x_1, Field6 sigma_x_2,                    // In
                                      Field7 sigma2_on_mu2_ip_1,                             // In
                                      Field8 sigma2_on_mu2_ip_2,                             // In
                                      bool l_const_Nc_in_cloud,                              // In
                                      Field9 mu_x_1_n, Field10 mu_x_2_n,                     // Out
                                      Field11 sigma_x_1_n, Field12 sigma_x_2_n               // Out
) {
  const auto hm_1_exec               = to_execution_space(hm_1);
  const auto hm_2_exec               = to_execution_space(hm_2);
  const auto Ncnm_exec               = to_execution_space(Ncnm);
  const auto mu_x_1_exec             = to_execution_space(mu_x_1);
  const auto mu_x_2_exec             = to_execution_space(mu_x_2);
  const auto sigma_x_1_exec          = to_execution_space(sigma_x_1);
  const auto sigma_x_2_exec          = to_execution_space(sigma_x_2);
  const auto sigma2_on_mu2_ip_1_exec = to_execution_space(sigma2_on_mu2_ip_1);
  const auto sigma2_on_mu2_ip_2_exec = to_execution_space(sigma2_on_mu2_ip_2);
  const auto mu_x_1_n_exec           = to_execution_space(mu_x_1_n, false);
  const auto mu_x_2_n_exec           = to_execution_space(mu_x_2_n, false);
  const auto sigma_x_1_n_exec        = to_execution_space(sigma_x_1_n, false);
  const auto sigma_x_2_n_exec        = to_execution_space(sigma_x_2_n, false);
  norm_transform_mean_stdev(nzt, ngrdcol, hydromet_dim, hm_metadata, hm_1_exec, hm_2_exec,
                            Ncnm_exec, pdf_dim, mu_x_1_exec, mu_x_2_exec, sigma_x_1_exec,
                            sigma_x_2_exec, sigma2_on_mu2_ip_1_exec, sigma2_on_mu2_ip_2_exec,
                            l_const_Nc_in_cloud, mu_x_1_n_exec, mu_x_2_n_exec, sigma_x_1_n_exec,
                            sigma_x_2_n_exec);
  from_execution_space(mu_x_1_n, mu_x_1_n_exec);
  from_execution_space(mu_x_2_n, mu_x_2_n_exec);
  from_execution_space(sigma_x_1_n, sigma_x_1_n_exec);
  from_execution_space(sigma_x_2_n, sigma_x_2_n_exec);
}
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7>
inline void denorm_transform_corr(int nzt, int ngrdcol, int pdf_dim,
                                  const hm_metadata_type& hm_metadata,                  // In
                                  Field0 sigma_x_1_n, Field1 sigma_x_2_n,               // In
                                  Field2 sigma2_on_mu2_ip_1, Field3 sigma2_on_mu2_ip_2, // In
                                  Field4 corr_array_1_n,                                // In
                                  Field5 corr_array_2_n,                                // In
                                  Field6 corr_array_1, Field7 corr_array_2              // Out
) {
  const auto sigma_x_1_n_exec        = to_execution_space(sigma_x_1_n);
  const auto sigma_x_2_n_exec        = to_execution_space(sigma_x_2_n);
  const auto sigma2_on_mu2_ip_1_exec = to_execution_space(sigma2_on_mu2_ip_1);
  const auto sigma2_on_mu2_ip_2_exec = to_execution_space(sigma2_on_mu2_ip_2);
  const auto corr_array_1_n_exec     = to_execution_space(corr_array_1_n);
  const auto corr_array_2_n_exec     = to_execution_space(corr_array_2_n);
  const auto corr_array_1_exec       = to_execution_space(corr_array_1, false);
  const auto corr_array_2_exec       = to_execution_space(corr_array_2, false);
  denorm_transform_corr(nzt, ngrdcol, pdf_dim, hm_metadata, sigma_x_1_n_exec, sigma_x_2_n_exec,
                        sigma2_on_mu2_ip_1_exec, sigma2_on_mu2_ip_2_exec, corr_array_1_n_exec,
                        corr_array_2_n_exec, corr_array_1_exec, corr_array_2_exec);
  from_execution_space(corr_array_1, corr_array_1_exec);
  from_execution_space(corr_array_2, corr_array_2_exec);
}
} // namespace clubb_kokkos
