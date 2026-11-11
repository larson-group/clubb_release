// Scalar KK numerical contracts; all implementations are C++.
#pragma once
#include "Microphys/microphys_host_types.hpp"
namespace clubb_kokkos::kk_engine {
inline double Dv_fnc(double order, double argument, const parameters_KK& kk);
inline real_t KK_Nrm_evap_upscaled_mean(
    real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1,
    real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n,
    real_t sigma_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n,
    real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_evap_coef,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t dt,
    const parameters_KK& kk);
inline real_t KK_Nrm_evap_local_mean(real_t KK_rrm_evap_tndcy, real_t Nrm, real_t rrm, real_t dt,
                                     const parameters_KK& kk);
inline real_t KK_Nrm_auto_mean(real_t KK_rrm_auto_tndcy, const parameters_KK& kk);
inline real_t KK_evap_local_mean(real_t mean_s, real_t rrm, real_t Nrm, real_t KK_evap_coef,
                                 const parameters_KK& kk);
inline real_t KK_auto_local_mean(real_t mean_s, real_t Ncm, real_t KK_auto_coef,
                                 const parameters_KK& kk);
inline real_t KK_accr_local_mean(real_t mean_s, real_t rrm, real_t KK_accr_coef,
                                 const parameters_KK& kk);
inline real_t KK_mvr_local_mean(real_t rrm, real_t Nrm, real_t KK_mvr_coef,
                                const parameters_KK& kk);
inline void KK_upscaled_means_driver(
    real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1,
    real_t mu_Nr_2, real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t mu_Ncn_1_n, real_t mu_Ncn_2_n, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_Ncn_1, real_t sigma_Ncn_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n,
    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t KK_evap_coef,
    real_t KK_auto_coef, real_t KK_accr_coef, real_t KK_mvr_coef, real_t& KK_evap_tndcy,
    real_t& KK_auto_tndcy, real_t& KK_accr_tndcy, real_t& KK_mean_vol_rad, const parameters_KK& kk);
inline real_t KK_evap_upscaled_mean(
    real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1,
    real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n,
    real_t sigma_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n,
    real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_evap_coef,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk);
inline real_t KK_auto_upscaled_mean(real_t mu_chi_1, real_t mu_chi_2, real_t mu_Ncn_1,
                                    real_t mu_Ncn_2, real_t mu_Ncn_1_n, real_t mu_Ncn_2_n,
                                    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                                    real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                                    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n,
                                    real_t KK_auto_coef, real_t mixt_frac, const parameters_KK& kk);
inline real_t KK_accr_upscaled_mean(real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
                                    real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                                    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_rr_1,
                                    real_t sigma_rr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
                                    real_t KK_accr_coef, real_t mixt_frac, real_t precip_frac_1,
                                    real_t precip_frac_2, const parameters_KK& kk);
inline real_t KK_mvr_upscaled_mean(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                   real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                   real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                   real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                   real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                   real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_mvr_coef,
                                   real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
                                   const parameters_KK& kk);
inline real_t trivar_NLL_mean_eq(real_t mu_chi_i, real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                 real_t mu_Nr_i_n, real_t sigma_chi_i, real_t sigma_rr_i,
                                 real_t sigma_Nr_i, real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                 real_t corr_chi_rr_i_n, real_t corr_chi_Nr_i_n,
                                 real_t corr_rr_Nr_i_n, real_t alpha_exp_in, real_t beta_exp_in,
                                 real_t gamma_exp_in, const parameters_KK& kk);
inline real_t bivar_NL_mean_eq(real_t mu_chi_i, real_t mu_y_i, real_t mu_y_i_n, real_t sigma_chi_i,
                               real_t sigma_y_i, real_t sigma_y_i_n, real_t corr_chi_y_i_n,
                               real_t y_tol, real_t alpha_exp_in, real_t beta_exp_in,
                               const parameters_KK& kk);
inline real_t bivar_LL_mean_eq(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n, real_t mu_Nr_i_n,
                               real_t sigma_rr_i, real_t sigma_Nr_i, real_t sigma_rr_i_n,
                               real_t sigma_Nr_i_n, real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                               real_t beta_exp_in, const parameters_KK& kk);
inline real_t variance_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_mean_vol_rad,
                              real_t KK_mvr_coef, real_t mixt_frac, real_t precip_frac_1,
                              real_t precip_frac_2, const parameters_KK& kk);
inline real_t G_T_p(real_t T_in_K, real_t p_in_Pa, int saturation_formula, const parameters_KK& kk);
inline double quadrivar_NNLL_covar(double mu_x1, double mu_x2, double mu_x3_n, double mu_x4_n,
                                   double sigma_x1, double sigma_x2, double sigma_x3_n,
                                   double sigma_x4_n, double rho_x1x2, double rho_x1x3_n,
                                   double rho_x1x4_n, double rho_x2x3_n, double rho_x2x4_n,
                                   double rho_x3x4_n, double x1_mean,
                                   double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                   double beta_exp, double gamma_exp, const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x1(double mu_x1, double mu_x2, double mu_x3_n,
                                            double mu_x4_n, double sigma_x2, double sigma_x3_n,
                                            double sigma_x4_n, double rho_x2x3_n, double rho_x2x4_n,
                                            double rho_x3x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x2(double mu_x1, double mu_x2, double mu_x3_n,
                                            double mu_x4_n, double sigma_x1, double sigma_x3_n,
                                            double sigma_x4_n, double rho_x1x3_n, double rho_x1x4_n,
                                            double rho_x3x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x3(double mu_x1, double mu_x2, double mu_x3,
                                            double mu_x4_n, double sigma_x1, double sigma_x2,
                                            double sigma_x4_n, double rho_x1x2, double rho_x1x4_n,
                                            double rho_x2x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                              double mu_x4_n, double sigma_x3_n, double sigma_x4_n,
                                              double rho_x3x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x1x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x2, double sigma_x4_n,
                                              double rho_x2x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x2x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x1, double sigma_x4_n,
                                              double rho_x1x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double sigma_x1, double sigma_x2,
                                              double rho_x1x2, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_cst_x1x2x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_cst_x1x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double sigma_x2, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_cst_x2x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk);
inline double quadrivar_NNLL_covar_const_all(double mu_x1, double mu_x2, double mu_x3, double mu_x4,
                                             double x1_mean, double x2_alpha_x3_beta_x4_gamma_mean,
                                             double alpha_exp, double beta_exp, double gamma_exp,
                                             const parameters_KK& kk);
inline double trivar_NNL_covar(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                               double sigma_x2, double sigma_x3_n, double rho_x1x2,
                               double rho_x1x3_n, double rho_x2x3_n, double x1_mean,
                               double x2_alpha_x3_beta_mean, double alpha_exp, double beta_exp,
                               const parameters_KK& kk);
inline double trivar_NNL_covar_const_x1(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x2,
                                        double sigma_x3_n, double rho_x2x3_n, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk);
inline double trivar_NNL_covar_const_x2(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                                        double sigma_x3_n, double rho_x1x3_n, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk);
inline double trivar_NNL_covar_const_x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x1,
                                        double sigma_x2, double rho_x1x2, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk);
inline double trivar_NNL_covar_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                          double sigma_x3_n, double x1_mean,
                                          double x2_alpha_x3_beta_mean, double alpha_exp,
                                          double beta_exp, const parameters_KK& kk);
inline double trivar_NNL_covar_const_x1x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x2,
                                          double x1_mean, double x2_alpha_x3_beta_mean,
                                          double alpha_exp, double beta_exp,
                                          const parameters_KK& kk);
inline double trivar_NNL_covar_const_x2x3(double mu_x1, double mu_x2, double mu_x3, double x1_mean,
                                          double x2_alpha_x3_beta_mean, double alpha_exp,
                                          double beta_exp, const parameters_KK& kk);
inline double trivar_NNL_covar_const_all(double mu_x1, double mu_x2, double mu_x3, double x1_mean,
                                         double x2_alpha_x3_beta_mean, double alpha_exp,
                                         double beta_exp, const parameters_KK& kk);
inline double trivar_NLL_mean(double mu_x1, double mu_x2_n, double mu_x3_n, double sigma_x1,
                              double sigma_x2_n, double sigma_x3_n, double rho_x1x2_n,
                              double rho_x1x3_n, double rho_x2x3_n, double alpha_exp,
                              double beta_exp, double gamma_exp, const parameters_KK& kk);
inline double trivar_NLL_mean_const_x1(double mu_x1, double mu_x2_n, double mu_x3_n,
                                       double sigma_x2_n, double sigma_x3_n, double rho_x2x3_n,
                                       double alpha_exp, double beta_exp, double gamma_exp,
                                       const parameters_KK& kk);
inline double trivar_NLL_mean_const_x2(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                                       double sigma_x3_n, double rho_x1x3_n, double alpha_exp,
                                       double beta_exp, double gamma_exp, const parameters_KK& kk);
inline double trivar_NLL_mean_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                         double sigma_x3_n, double alpha_exp, double beta_exp,
                                         double gamma_exp, const parameters_KK& kk);
inline double trivar_NLL_mean_const_x2x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x1,
                                         double alpha_exp, double beta_exp, double gamma_exp,
                                         const parameters_KK& kk);
inline double trivar_NLL_mean_const_all(double mu_x1, double mu_x2, double mu_x3, double alpha_exp,
                                        double beta_exp, double gamma_exp, const parameters_KK& kk);
inline double bivar_NL_mean(double mu_x1, double mu_x2_n, double sigma_x1, double sigma_x2_n,
                            double rho_x1x2_n, double alpha_exp, double beta_exp,
                            const parameters_KK& kk);
inline double bivar_NL_mean_const_x1(double mu_x1, double mu_x2_n, double sigma_x2_n,
                                     double alpha_exp, double beta_exp, const parameters_KK& kk);
inline double bivar_NL_mean_const_x2(double mu_x1, double mu_x2, double sigma_x1, double alpha_exp,
                                     double beta_exp, const parameters_KK& kk);
inline double bivar_NL_mean_const_all(double mu_x1, double mu_x2, double alpha_exp, double beta_exp,
                                      const parameters_KK& kk);
inline double bivar_LL_mean(double mu_x1_n, double mu_x2_n, double sigma_x1_n, double sigma_x2_n,
                            double rho_x1x2_n, double alpha_exp, double beta_exp,
                            const parameters_KK& kk);
inline double bivar_LL_mean_const_x1(double mu_x1, double mu_x2_n, double sigma_x2_n,
                                     double alpha_exp, double beta_exp, const parameters_KK& kk);
inline double bivar_LL_mean_const_all(double mu_x1, double mu_x2, double alpha_exp, double beta_exp,
                                      const parameters_KK& kk);
inline void KK_upscaled_covar_driver(
    real_t w_mean, real_t rtm, real_t thlm, real_t exner, real_t rrainm, real_t Nrm, real_t mu_w_1,
    real_t mu_w_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_eta_1, real_t mu_eta_2,
    real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_Ncn_1,
    real_t mu_Ncn_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
    real_t mu_Ncn_1_n, real_t mu_Ncn_2_n, real_t sigma_w_1, real_t sigma_w_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_rr_1,
    real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_Ncn_1, real_t sigma_Ncn_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n, real_t corr_w_chi_1, real_t corr_w_chi_2,
    real_t corr_w_rr_1_n, real_t corr_w_rr_2_n, real_t corr_w_Nr_1_n, real_t corr_w_Nr_2_n,
    real_t corr_w_Ncn_1_n, real_t corr_w_Ncn_2_n, real_t corr_chi_eta_1, real_t corr_chi_eta_2,
    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n,
    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t corr_eta_rr_1_n,
    real_t corr_eta_rr_2_n, real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_eta_Ncn_1_n,
    real_t corr_eta_Ncn_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t mixt_frac,
    real_t precip_frac_1, real_t precip_frac_2, real_t KK_evap_coef, real_t KK_auto_coef,
    real_t KK_accr_coef, real_t KK_evap_tndcy, real_t KK_auto_tndcy, real_t KK_accr_tndcy,
    real_t mu_rt_1, real_t mu_rt_2, real_t mu_thl_1, real_t mu_thl_2, real_t crt1, real_t crt2,
    real_t cthl1, real_t cthl2, real_t& wprtp_mc_src_tndcy, real_t& wpthlp_mc_src_tndcy,
    real_t& rtp2_mc_src_tndcy, real_t& thlp2_mc_src_tndcy, real_t& rtpthlp_mc_src_tndcy,
    real_t& w_KK_evap_covar, real_t& rt_KK_evap_covar, real_t& thl_KK_evap_covar,
    real_t& w_KK_auto_covar, real_t& rt_KK_auto_covar, real_t& thl_KK_auto_covar,
    real_t& w_KK_accr_covar, real_t& rt_KK_accr_covar, real_t& thl_KK_accr_covar,
    const parameters_KK& kk);
inline real_t covar_x_KK_evap(
    real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2,
    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
    real_t mu_Nr_2_n, real_t sigma_x_1, real_t sigma_x_2, real_t sigma_chi_1, real_t sigma_chi_2,
    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_x_chi_1,
    real_t corr_x_chi_2, real_t corr_x_rr_1_n, real_t corr_x_rr_2_n, real_t corr_x_Nr_1_n,
    real_t corr_x_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n,
    real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t x_mean,
    real_t KK_evap_tndcy, real_t KK_evap_coef, real_t x_tol, real_t mixt_frac, real_t precip_frac_1,
    real_t precip_frac_2, const parameters_KK& kk);
inline real_t covar_rt_KK_evap(
    real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
    real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n, real_t corr_eta_rr_2_n,
    real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
    real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t rtm, real_t mu_rt_1,
    real_t mu_rt_2, real_t KK_evap_tndcy, real_t KK_evap_coef, real_t eta_tol, real_t crt1,
    real_t crt2, const parameters_KK& kk);
inline real_t covar_thl_KK_evap(
    real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
    real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n, real_t corr_eta_rr_2_n,
    real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
    real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t thlm, real_t mu_thl_1,
    real_t mu_thl_2, real_t KK_evap_tndcy, real_t KK_evap_coef, real_t eta_tol, real_t cthl1,
    real_t cthl2, const parameters_KK& kk);
inline real_t covar_x_KK_auto(real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2,
                              real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                              real_t mu_Ncn_2_n, real_t sigma_x_1, real_t sigma_x_2,
                              real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                              real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                              real_t corr_x_chi_1, real_t corr_x_chi_2, real_t corr_x_Ncn_1_n,
                              real_t corr_x_Ncn_2_n, real_t corr_chi_Ncn_1_n,
                              real_t corr_chi_Ncn_2_n, real_t x_mean, real_t KK_auto_tndcy,
                              real_t KK_auto_coef, real_t x_tol, real_t mixt_frac,
                              const parameters_KK& kk);
inline real_t covar_rt_KK_auto(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                               real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                               real_t mu_Ncn_2_n, real_t sigma_eta_1, real_t sigma_eta_2,
                               real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                               real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                               real_t corr_eta_chi_1, real_t corr_eta_chi_2,
                               real_t corr_eta_Ncn_1_n, real_t corr_eta_Ncn_2_n,
                               real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t rtm,
                               real_t mu_rt_1, real_t mu_rt_2, real_t KK_auto_tndcy,
                               real_t KK_auto_coef, real_t eta_tol, real_t crt1, real_t crt2,
                               real_t mixt_frac, const parameters_KK& kk);
inline real_t covar_thl_KK_auto(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                                real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                                real_t mu_Ncn_2_n, real_t sigma_eta_1, real_t sigma_eta_2,
                                real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                                real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                                real_t corr_eta_chi_1, real_t corr_eta_chi_2,
                                real_t corr_eta_Ncn_1_n, real_t corr_eta_Ncn_2_n,
                                real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t thlm,
                                real_t mu_thl_1, real_t mu_thl_2, real_t KK_auto_tndcy,
                                real_t KK_auto_coef, real_t eta_tol, real_t cthl1, real_t cthl2,
                                real_t mixt_frac, const parameters_KK& kk);
inline real_t covar_x_KK_accr(real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2,
                              real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                              real_t sigma_x_1, real_t sigma_x_2, real_t sigma_chi_1,
                              real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t corr_x_chi_1,
                              real_t corr_x_chi_2, real_t corr_x_rr_1_n, real_t corr_x_rr_2_n,
                              real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t x_mean,
                              real_t KK_accr_tndcy, real_t KK_accr_coef, real_t x_tol,
                              real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
                              const parameters_KK& kk);
inline real_t covar_rt_KK_accr(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                               real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                               real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
                               real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2,
                               real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t corr_eta_chi_1,
                               real_t corr_eta_chi_2, real_t corr_eta_rr_1_n,
                               real_t corr_eta_rr_2_n, real_t corr_chi_rr_1_n,
                               real_t corr_chi_rr_2_n, real_t rtm, real_t mu_rt_1, real_t mu_rt_2,
                               real_t KK_accr_tndcy, real_t KK_accr_coef, real_t eta_tol,
                               real_t crt1, real_t crt2, real_t mixt_frac, real_t precip_frac_1,
                               real_t precip_frac_2, const parameters_KK& kk);
inline real_t
covar_thl_KK_accr(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                  real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                  real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1, real_t sigma_chi_2,
                  real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                  real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n,
                  real_t corr_eta_rr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
                  real_t thlm, real_t mu_thl_1, real_t mu_thl_2, real_t KK_accr_tndcy,
                  real_t KK_accr_coef, real_t eta_tol, real_t cthl1, real_t cthl2, real_t mixt_frac,
                  real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk);
inline real_t quadrivar_NNLL_covar_eq(
    real_t mu_x_i, real_t mu_chi_i, real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
    real_t mu_Nr_i_n, real_t sigma_x_i, real_t sigma_chi_i, real_t sigma_rr_i, real_t sigma_Nr_i,
    real_t sigma_rr_i_n, real_t sigma_Nr_i_n, real_t corr_x_chi_i, real_t corr_x_rr_i_n,
    real_t corr_x_Nr_i_n, real_t corr_chi_rr_i_n, real_t corr_chi_Nr_i_n, real_t corr_rr_Nr_i_n,
    real_t x_mean, real_t mc_tndcy_mean, real_t mc_coef, real_t x_tol, real_t alpha_exp_in,
    real_t beta_exp_in, real_t gamma_exp_in, const parameters_KK& kk);
inline real_t trivar_NNL_covar_eq(real_t mu_x_i, real_t mu_chi_i, real_t mu_y_i, real_t mu_y_i_n,
                                  real_t sigma_x_i, real_t sigma_chi_i, real_t sigma_y_i,
                                  real_t sigma_y_i_n, real_t corr_x_chi_i, real_t corr_x_y_i_n,
                                  real_t corr_chi_y_i_n, real_t x_mean, real_t mc_tndcy_mean,
                                  real_t mc_coef, real_t x_tol, real_t y_tol, real_t alpha_exp_in,
                                  real_t beta_exp_in, const parameters_KK& kk);
inline void KK_sed_vel_covars(real_t rrm, real_t rr_1, real_t rr_2, real_t Nrm, real_t Nr_1,
                              real_t Nr_2, real_t KK_mean_vol_rad, real_t mu_rr_1, real_t mu_rr_2,
                              real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                              real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_rr_1,
                              real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
                              real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n,
                              real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                              real_t KK_mvr_coef, real_t mixt_frac, real_t precip_frac_1,
                              real_t precip_frac_2, real_t& Vrrprrp_impc, real_t& Vrrprrp_expc,
                              real_t& VNrpNrp_impc, real_t& VNrpNrp_expc, real_t& rr_KK_mvr_covar,
                              real_t& Nr_KK_mvr_covar, const parameters_KK& kk);
inline real_t covar_rr_KK_mvr_coefA(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                    real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                    real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                    real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                    real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                                    real_t KK_mean_vol_rad, real_t KK_mvr_coef,
                                    const parameters_KK& kk);
inline real_t covar_rr_KK_mvr_termB(real_t rr_1, real_t rr_2, real_t mu_rr_1, real_t mu_rr_2,
                                    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n,
                                    real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
                                    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
                                    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n,
                                    real_t corr_rr_Nr_2_n, real_t KK_mvr_coef, real_t mixt_frac,
                                    const parameters_KK& kk);
inline real_t covar_rr_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t rrm,
                              real_t KK_mean_vol_rad, real_t KK_mvr_coef, real_t mixt_frac,
                              real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk);
inline real_t covar_Nr_KK_mvr_coefA(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                    real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                    real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                    real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                    real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                                    real_t KK_mean_vol_rad, real_t KK_mvr_coef,
                                    const parameters_KK& kk);
inline real_t covar_Nr_KK_mvr_termB(real_t Nr_1, real_t Nr_2, real_t mu_rr_1, real_t mu_rr_2,
                                    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n,
                                    real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
                                    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
                                    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n,
                                    real_t corr_rr_Nr_2_n, real_t KK_mvr_coef, real_t mixt_frac,
                                    const parameters_KK& kk);
inline real_t covar_Nr_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t Nrm,
                              real_t KK_mean_vol_rad, real_t KK_mvr_coef, real_t mixt_frac,
                              real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk);
inline real_t bivar_LL_covar_partial_rr(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                        real_t mu_Nr_i_n, real_t sigma_rr_i, real_t sigma_Nr_i,
                                        real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                        real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                                        real_t beta_exp_in, const parameters_KK& kk);
inline real_t bivar_LL_covar_partial_Nr(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                        real_t mu_Nr_i_n, real_t sigma_rr_i, real_t sigma_Nr_i,
                                        real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                        real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                                        real_t beta_exp_in, const parameters_KK& kk);
inline double bivar_LL_covar_partial(double mu_x1_n, double mu_x2_n, double sigma_x1_n,
                                     double sigma_x2_n, double rho_x1x2_n, double alpha_exp,
                                     double beta_exp, const parameters_KK& kk);
inline double bivar_LL_covar_const_x2_partial(double mu_x1_n, double mu_x2, double sigma_x1_n,
                                              double alpha_exp, double beta_exp,
                                              const parameters_KK& kk);
} // namespace clubb_kokkos::kk_engine
