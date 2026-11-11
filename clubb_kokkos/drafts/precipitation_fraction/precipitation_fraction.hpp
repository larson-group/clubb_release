// Diagnose overall and component precipitation fractions for CLUBB's PDF.
// precip_fraction retains warm/ice options, both physical grid directions,
// species-ordered boosts, Native stats and full first-failure diagnostics.
// Fields use column-first LayoutLeft; vertical maximum recurrences stay serial
// within each independent column. Scratch is per call; no mutable global cache.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
inline constexpr int precip_frac_calc_type = 2; // Native public parameter.
// gr: Grid type variable
// nzt: Number of model thermodynamic vertical grid levels
// ngrdcol: Number of grid columns
// hydromet_dim: Number of hydrometeor species
// hydromet: Mean of hydrometeor, hm (overall)           [units vary]
// cloud_frac: Cloud fraction (overall)                      [-]
// cloud_frac_1: Cloud fraction (1st PDF component)            [-]
// l_mix_rat_hm: if true, then the quantity is a hydrometeor mixing ratio
// l_frozen_hm: if true, then the hydrometeor is frozen; otherwise liquid
// hydromet_tol: Tolerance values for all hydrometeors    [units vary]
// cloud_frac_2: Cloud fraction (2nd PDF component)            [-]
// ice_supersat_frac: Ice supersaturation fraction (overall)        [-]
// ice_supersat_frac_1: Ice supersaturation fraction (1st PDF comp.)  [-]
// ice_supersat_frac_2: Ice supersaturation fraction (2nd PDF comp.)  [-]
// mixt_frac: Mixture fraction                              [-]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// err_info: err_info struct containing err_code and err_header
// precip_frac: Precipitation fraction (overall)               [-]
// precip_frac_1: Precipitation fraction (1st PDF component)     [-]
// precip_frac_2: Precipitation fraction (2nd PDF component)     [-]
// precip_frac_tol: Minimum precip. frac. when hydromet. are present  [-]
// stats: Native statistics collector.
void precip_fraction(const grid& gr, int nzt, int ngrdcol, int hydromet_dim,       // In
                     Real3D hydromet, Real2D cloud_frac, Real2D cloud_frac_1,      // In
                     Bool1D l_mix_rat_hm, Bool1D l_frozen_hm, Real1D hydromet_tol, // In
                     Real2D cloud_frac_2, Real2D ice_supersat_frac,                // In
                     Real2D ice_supersat_frac_1, Real2D ice_supersat_frac_2,       // In
                     Real2D mixt_frac, Real2D clubb_params,                        // In
                     err_info_type& err_info,                                      // InOut
                     Real2D precip_frac,                                           // Out
                     Real2D precip_frac_1,                                         // Out
                     Real2D precip_frac_2,                                         // Out
                     Real1D precip_frac_tol,                                       // Out
                     f2k::stats_type& stats                                        // InOut
);
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15>
inline void precip_fraction(const grid& gr, int nzt, int ngrdcol, int hydromet_dim,       // In
                            Field0 hydromet, Field1 cloud_frac, Field2 cloud_frac_1,      // In
                            Field3 l_mix_rat_hm, Field4 l_frozen_hm, Field5 hydromet_tol, // In
                            Field6 cloud_frac_2, Field7 ice_supersat_frac,                // In
                            Field8 ice_supersat_frac_1, Field9 ice_supersat_frac_2,       // In
                            Field10 mixt_frac, Field11 clubb_params,                      // In
                            err_info_type& err_info,                                      // InOut
                            Field12 precip_frac,                                          // Out
                            Field13 precip_frac_1,                                        // Out
                            Field14 precip_frac_2,                                        // Out
                            Field15 precip_frac_tol,                                      // Out
                            f2k::stats_type& stats                                        // InOut
) {
  const auto hydromet_exec = to_execution_space(hydromet);
  const auto cloud_frac_exec = to_execution_space(cloud_frac);
  const auto cloud_frac_1_exec = to_execution_space(cloud_frac_1);
  const auto l_mix_rat_hm_exec = to_execution_space(l_mix_rat_hm);
  const auto l_frozen_hm_exec = to_execution_space(l_frozen_hm);
  const auto hydromet_tol_exec = to_execution_space(hydromet_tol);
  const auto cloud_frac_2_exec = to_execution_space(cloud_frac_2);
  const auto ice_supersat_frac_exec = to_execution_space(ice_supersat_frac);
  const auto ice_supersat_frac_1_exec = to_execution_space(ice_supersat_frac_1);
  const auto ice_supersat_frac_2_exec = to_execution_space(ice_supersat_frac_2);
  const auto mixt_frac_exec = to_execution_space(mixt_frac);
  const auto clubb_params_exec = to_execution_space(clubb_params);
  const auto precip_frac_exec = to_execution_space(precip_frac, false);
  const auto precip_frac_1_exec = to_execution_space(precip_frac_1, false);
  const auto precip_frac_2_exec = to_execution_space(precip_frac_2, false);
  const auto precip_frac_tol_exec = to_execution_space(precip_frac_tol, false);
  precip_fraction(gr, nzt, ngrdcol, hydromet_dim, hydromet_exec, cloud_frac_exec, cloud_frac_1_exec,
                  l_mix_rat_hm_exec, l_frozen_hm_exec, hydromet_tol_exec, cloud_frac_2_exec,
                  ice_supersat_frac_exec, ice_supersat_frac_1_exec, ice_supersat_frac_2_exec,
                  mixt_frac_exec, clubb_params_exec, err_info, precip_frac_exec, precip_frac_1_exec,
                  precip_frac_2_exec, precip_frac_tol_exec, stats);
  from_execution_space(precip_frac, precip_frac_exec);
  from_execution_space(precip_frac_1, precip_frac_1_exec);
  from_execution_space(precip_frac_2, precip_frac_2_exec);
  from_execution_space(precip_frac_tol, precip_frac_tol_exec);
}
} // namespace clubb_kokkos
