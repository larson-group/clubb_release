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
// src/CLUBB_core/diffusion.F90:diffusion_zm_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void diffusion_zm_lhs(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                      const Array0& K_zt, const Array1& K_zm, const Array2& nu,
                      const Array3& invrs_rho_ds_zm, const Array4& rho_ds_zt, const Array5& lhs);
// src/CLUBB_core/mean_adv.F90:term_ma_zm_lhs; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void term_ma_zm_lhs(int nzm, int nzt, int ngrdcol, const Array0& wm_zm, const Array1& invrs_dzm,
                    const Array2& weights_zm2zt, const Array3& lhs_ma);
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
// src/CLUBB_core/turbulent_adv_pdf.F90:xpyp_term_ta_pdf_rhs_godunov; all arguments required.
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void xpyp_term_ta_pdf_rhs_godunov(int nzm, int nzt, int ngrdcol, const clubb_kokkos::grid& gr,
                                  const Array0& term_wpxpyp_explicit_zm,
                                  const Array1& invrs_rho_ds_zm, const Array2& sgn_turbulent_vel,
                                  const Array3& rho_ds_zm, const Array4& rhs_ta);
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
// src/CLUBB_core/clip_explicit.F90:clip_variance; Native optional threshold_hi retains presence
// semantics.
template <class Array0, class Array1>
void clip_variance(int nzm, int ngrdcol, const clubb_kokkos::grid& gr, int solve_type,
                   clubb_kokkos::real_t dt, const Array0& threshold_lo, stats_type& stats,
                   const Array1& xp2, const clubb_kokkos::real_t* threshold_hi = nullptr);
// src/CLUBB_core/stats_netcdf.F90:stats_update_budget_2d; all arguments required.
template <class Array0>
void stats_update_budget(const std::string& name, const Array0& values, stats_type& stats);
// Native stats_finalize_budget_2d optional l_count_sample; explicit-value overload.
template <class Array0>
void stats_finalize_budget(const std::string& name, const Array0& values, stats_type& stats,
                           bool l_count_sample);
} // namespace f2k
