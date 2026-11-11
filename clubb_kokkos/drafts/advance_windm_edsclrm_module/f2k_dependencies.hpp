// Declaration-only Native lower boundaries for this isolated draft.
// Root must add missing definitions through the existing f2k generator owner,
// then remove this staging header. These declarations are not bridge implementations.
#pragma once
#include "f2k/procedure_bridge.hpp"
namespace f2k {
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
// src/CLUBB_core/clip_explicit.F90:clip_covar; all arguments required.
template <class Array0, class Array1, class Array2, class Array3>
void clip_covar(int nzm, int ngrdcol, int solve_type, const Array0& xp2, const Array1& yp2,
                const Array2& xpyp, const Array3& xpyp_chnge);
// src/CLUBB_core/fill_holes.F90:fill_holes_vertical_api; all arguments required.
template <class Array0, class Array1, class Array2>
void fill_holes_vertical_api(int nz, int ngrdcol, clubb_kokkos::real_t threshold,
                             int lower_hf_level, int upper_hf_level, const Array0& dz,
                             const Array1& rho_ds, int grid_dir_indx, int fill_holes_type,
                             const Array2& field);
// src/CLUBB_core/stats_netcdf.F90:stats_update_budget_2d; all arguments required.
template <class Array0>
void stats_update_budget(const std::string& name, const Array0& values, stats_type& stats);
} // namespace f2k
