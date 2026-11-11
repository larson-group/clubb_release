// Declaration-only Native lower boundaries for this isolated draft.
// Root must add missing definitions through the existing f2k generator owner,
// then remove this staging header. These declarations are not bridge implementations.
#pragma once
#include "f2k/procedure_bridge.hpp"
namespace f2k {
// src/CLUBB_core/advance_helper_module.F90:calc_wp3_on_wp2; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void calc_wp3_on_wp2(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr, const Array0& wp2,
                     const Array1& wp3, const Array2& wp3_on_wp2, const Array3& wp3_on_wp2_zt);
// src/CLUBB_core/mono_flux_limiter.F90:calc_turb_adv_range; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void calc_turb_adv_range(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                         clubb_kokkos::real_t dt, const Array0& w_1_zm, const Array1& w_2_zm,
                         const Array2& varnce_w_1_zm, const Array3& varnce_w_2_zm,
                         const Array4& mixt_frac_zm, stats_type& stats,
                         const Array5& low_lev_effect, const Array6& high_lev_effect);
// src/CLUBB_core/diffusion.F90:diffusion_zm_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void diffusion_zm_lhs(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                      const Array0& K_zt, const Array1& K_zm, const Array2& nu,
                      const Array3& invrs_rho_ds_zm, const Array4& rho_ds_zt, const Array5& lhs);
// src/CLUBB_core/mean_adv.F90:term_ma_zm_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void term_ma_zm_lhs(int nzm, int nzt, int ngrdcol, const Array0& wm_zm, const Array1& invrs_dzm,
                    const Array2& weights_zm2zt, const Array3& lhs_ma);
// src/CLUBB_core/diffusion.F90:diffusion_zt_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void diffusion_zt_lhs(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                      const Array0& K_zm, const Array1& K_zt, const Array2& nu,
                      const Array3& invrs_rho_ds_ztzxt, const Array4& rho_ds_zm, const Array5& lhs);
// src/CLUBB_core/mean_adv.F90:term_ma_zt_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void term_ma_zt_lhs(int nzm, int nzt, int ngrdcol, const Array0& wm_zt, const Array1& weights_zt2zm,
                    const Array2& invrs_dzt, const Array3& invrs_dzm, bool l_upwind_xm_ma,
                    clubb_kokkos::real_t grid_dir, const Array4& lhs_ma);
// src/CLUBB_core/turbulent_adv_pdf.F90:xpyp_term_ta_pdf_rhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void xpyp_term_ta_pdf_rhs(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                          const Array0& term_wpxpyp_explicit, const Array1& rho_ds_zt,
                          const Array2& rho_ds_zm, const Array3& invrs_rho_ds_zm,
                          bool l_upwind_xpyp_turbulent_adv, const Array4& sgn_turbulent_vel,
                          const Array5& term_wpxpyp_explicit_zm, const Array6& rhs_ta);
// src/CLUBB_core/turbulent_adv_pdf.F90:xpyp_term_ta_pdf_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void xpyp_term_ta_pdf_lhs(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                          const Array0& coef_wpxpyp_implicit, const Array1& rho_ds_zt,
                          const Array2& rho_ds_zm, const Array3& invrs_rho_ds_zm,
                          bool l_upwind_xpyp_turbulent_adv, const Array4& sgn_turbulent_vel,
                          const Array5& coef_wpxpyp_implicit_zm, const Array6& lhs_ta);
// src/CLUBB_core/turbulent_adv_pdf.F90:xpyp_term_ta_pdf_lhs_godunov; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void xpyp_term_ta_pdf_lhs_godunov(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                                  const Array0& coef_wpxpyp_implicit, const Array1& invrs_rho_ds_zm,
                                  const Array2& rho_ds_zm, const Array3& lhs_ta);
// src/CLUBB_core/mono_flux_limiter.F90:monotonic_turbulent_flux_limit; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11>
void monotonic_turbulent_flux_limit(
    int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr, int solve_type,
    clubb_kokkos::real_t dt, const Array0& xm_old, const Array1& xp2, const Array2& wm_zt,
    const Array3& xm_forcing, const Array4& rho_ds_zm, const Array5& rho_ds_zt,
    const Array6& invrs_rho_ds_zm, const Array7& invrs_rho_ds_zt,
    clubb_kokkos::real_t xp2_threshold, clubb_kokkos::real_t xm_tol, bool l_implemented,
    const Array8& low_lev_effect, const Array9& high_lev_effect, int tridiag_solve_method,
    bool l_upwind_xm_ma, bool l_mono_flux_lim_spikefix, stats_type& stats, const Array10& xm,
    const Array11& wpxp, clubb_kokkos::err_info_type& err_info);
// src/CLUBB_core/pos_definite_module.F90:pos_definite_adj; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void pos_definite_adj(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                      clubb_kokkos::real_t dt, const Array0& field_np1, const Array1& flux_np1,
                      const Array2& field_n, const Array3& field_pd, const Array4& flux_pd);
// src/CLUBB_core/fill_holes.F90:fill_holes_vertical_api; all arguments required.
template <class Array0, class Array1, class Array2>
void fill_holes_vertical_api(int nz, int ngrdcol, clubb_kokkos::real_t threshold,
                             int lower_hf_level, int upper_hf_level, const Array0& dz,
                             const Array1& rho_ds, int grid_dir_indx, int fill_holes_type,
                             const Array2& field);
// src/CLUBB_core/clip_explicit.F90:clip_covar; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void clip_covar(int nzm, int ngrdcol, int solve_type, const Array0& xp2, const Array1& yp2,
                const Array2& xpyp, const Array3& xpyp_chnge);
// src/CLUBB_core/stats_netcdf.F90:stats_update_budget_2d; all arguments required.
template <class Array0>
void stats_update_budget(const std::string& name, const Array0& values, stats_type& stats);
} // namespace f2k
