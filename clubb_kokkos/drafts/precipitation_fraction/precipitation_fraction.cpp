// Precipitation fraction port of src/CLUBB_core/precipitation_fraction.F90.
#include "precipitation_fraction.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
namespace clubb_kokkos {
namespace {
void component_precip_frac_weighted(const grid& gr, int nzt, int ngrdcol, int hydromet_dim, // In
                                    Bool1D l_frozen_hm, Real1D hydromet_tol,                // In
                                    Real3D hydromet, Real2D precip_frac,                    // In
                                    Real2D cloud_frac_1, Real2D cloud_frac_2,               // In
                                    Real2D ice_supersat_frac_1,                             // In
                                    Real2D ice_supersat_frac_2, Real2D mixt_frac,           // In
                                    Real1D precip_frac_tol,                                 // In
                                    Real2D precip_frac_1, Real2D precip_frac_2              // Out
);
void component_precip_frac_specify(int nzt, int ngrdcol, int hydromet_dim,
                                   Real1D hydromet_tol,                       // In
                                   Real1D upsilon_precip_frac_rat,            // In
                                   Real3D hydromet, Real2D precip_frac,       // In
                                   Real2D mixt_frac, Real1D precip_frac_tol,  // In
                                   Real2D precip_frac_1, Real2D precip_frac_2 // Out
);
void precip_frac_assert_check(int nzt, int hydromet_dim, Real1D hydromet_tol,        // In
                              Real2D hydromet, Real1D mixt_frac, Real1D precip_frac, // In
                              Real1D precip_frac_1, Real1D precip_frac_2,            // In
                              real_t precip_frac_tol,                                // In
                              err_info_type& err_info                                // InOut
);
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11>
inline void component_precip_frac_weighted(const grid& gr, int nzt, int ngrdcol,
                                           int hydromet_dim,                             // In
                                           Field0 l_frozen_hm, Field1 hydromet_tol,      // In
                                           Field2 hydromet, Field3 precip_frac,          // In
                                           Field4 cloud_frac_1, Field5 cloud_frac_2,     // In
                                           Field6 ice_supersat_frac_1,                   // In
                                           Field7 ice_supersat_frac_2, Field8 mixt_frac, // In
                                           Field9 precip_frac_tol,                       // In
                                           Field10 precip_frac_1, Field11 precip_frac_2  // Out
) {
  const auto l_frozen_hm_exec = to_execution_space(l_frozen_hm);
  const auto hydromet_tol_exec = to_execution_space(hydromet_tol);
  const auto hydromet_exec = to_execution_space(hydromet);
  const auto precip_frac_exec = to_execution_space(precip_frac);
  const auto cloud_frac_1_exec = to_execution_space(cloud_frac_1);
  const auto cloud_frac_2_exec = to_execution_space(cloud_frac_2);
  const auto ice_supersat_frac_1_exec = to_execution_space(ice_supersat_frac_1);
  const auto ice_supersat_frac_2_exec = to_execution_space(ice_supersat_frac_2);
  const auto mixt_frac_exec = to_execution_space(mixt_frac);
  const auto precip_frac_tol_exec = to_execution_space(precip_frac_tol);
  const auto precip_frac_1_exec = to_execution_space(precip_frac_1, false);
  const auto precip_frac_2_exec = to_execution_space(precip_frac_2, false);
  component_precip_frac_weighted(gr, nzt, ngrdcol, hydromet_dim, l_frozen_hm_exec,
                                 hydromet_tol_exec, hydromet_exec, precip_frac_exec,
                                 cloud_frac_1_exec, cloud_frac_2_exec, ice_supersat_frac_1_exec,
                                 ice_supersat_frac_2_exec, mixt_frac_exec, precip_frac_tol_exec,
                                 precip_frac_1_exec, precip_frac_2_exec);
  from_execution_space(precip_frac_1, precip_frac_1_exec);
  from_execution_space(precip_frac_2, precip_frac_2_exec);
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7>
inline void component_precip_frac_specify(int nzt, int ngrdcol, int hydromet_dim,
                                          Field0 hydromet_tol,                       // In
                                          Field1 upsilon_precip_frac_rat,            // In
                                          Field2 hydromet, Field3 precip_frac,       // In
                                          Field4 mixt_frac, Field5 precip_frac_tol,  // In
                                          Field6 precip_frac_1, Field7 precip_frac_2 // Out
) {
  const auto hydromet_tol_exec = to_execution_space(hydromet_tol);
  const auto upsilon_precip_frac_rat_exec = to_execution_space(upsilon_precip_frac_rat);
  const auto hydromet_exec = to_execution_space(hydromet);
  const auto precip_frac_exec = to_execution_space(precip_frac);
  const auto mixt_frac_exec = to_execution_space(mixt_frac);
  const auto precip_frac_tol_exec = to_execution_space(precip_frac_tol);
  const auto precip_frac_1_exec = to_execution_space(precip_frac_1, false);
  const auto precip_frac_2_exec = to_execution_space(precip_frac_2, false);
  component_precip_frac_specify(nzt, ngrdcol, hydromet_dim, hydromet_tol_exec,
                                upsilon_precip_frac_rat_exec, hydromet_exec, precip_frac_exec,
                                mixt_frac_exec, precip_frac_tol_exec, precip_frac_1_exec,
                                precip_frac_2_exec);
  from_execution_space(precip_frac_1, precip_frac_1_exec);
  from_execution_space(precip_frac_2, precip_frac_2_exec);
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5>
inline void precip_frac_assert_check(int nzt, int hydromet_dim, Field0 hydromet_tol,        // In
                                     Field1 hydromet, Field2 mixt_frac, Field3 precip_frac, // In
                                     Field4 precip_frac_1, Field5 precip_frac_2,            // In
                                     real_t precip_frac_tol,                                // In
                                     err_info_type& err_info                                // InOut
) {
  const auto hydromet_tol_exec = to_execution_space(hydromet_tol);
  const auto hydromet_exec = to_execution_space(hydromet);
  const auto mixt_frac_exec = to_execution_space(mixt_frac);
  const auto precip_frac_exec = to_execution_space(precip_frac);
  const auto precip_frac_1_exec = to_execution_space(precip_frac_1);
  const auto precip_frac_2_exec = to_execution_space(precip_frac_2);
  precip_frac_assert_check(nzt, hydromet_dim, hydromet_tol_exec, hydromet_exec, mixt_frac_exec,
                           precip_frac_exec, precip_frac_1_exec, precip_frac_2_exec,
                           precip_frac_tol, err_info);
}

} // namespace
// =============================================================================
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
) {
  // Description:
  // Determines (overall) precipitation fraction over the horizontal domain, as
  // well as the precipitation fraction within each PDF component, at every
  // vertical grid level.
  // References:
  //-----------------------------------------------------------------------

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

  constexpr real_t max_hm_ip_comp_mean = real_t(0.0025); // [kg/kg]
  constexpr real_t precip_frac_tol_coef = real_t(0.1);   // Coefficient for precip_frac_tol
  int j{};                                               // Loop indices
  int ivar{};
  // Reduce immutable species flags once per call; no floating reduction.
  int frozen_hm_count = 0;
  Kokkos::parallel_reduce(
      "precip_fraction::frozen_species", hydromet_dim,
      KOKKOS_LAMBDA(const int species, int& count) { count += l_frozen_hm(species) != 0; },
      frozen_hm_count);
  const bool frozen_hm_present = frozen_hm_count != 0;
  const auto l_mix_rat_hm_host = to_host_space(l_mix_rat_hm);
  // Capture component Views/scalars only; owning UDTs remain on the host.
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zt = gr.k_lb_zt;

  // Initialize the precipitation fraction variables (precip_frac,
  // precip_frac_1, and precip_frac_2) to 0.
  {
    const int slice_extent_0 = precip_frac.extent_int(0);
    const int slice_extent_1 = precip_frac.extent_int(1);
    Kokkos::parallel_for(
        "precip_fraction::array_1", slice_extent_0 * slice_extent_1, KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = cell / (1) % slice_extent_0 + 1;
          const int slice_1 = cell / (1 * slice_extent_0) % slice_extent_1 + 1;
          precip_frac(slice_0 - 1, slice_1 - 1) = zero;
        });
  }
  {
    const int slice_extent_0 = precip_frac_1.extent_int(0);
    const int slice_extent_1 = precip_frac_1.extent_int(1);
    Kokkos::parallel_for(
        "precip_fraction::array_2", slice_extent_0 * slice_extent_1, KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = cell / (1) % slice_extent_0 + 1;
          const int slice_1 = cell / (1 * slice_extent_0) % slice_extent_1 + 1;
          precip_frac_1(slice_0 - 1, slice_1 - 1) = zero;
        });
  }
  {
    const int slice_extent_0 = precip_frac_2.extent_int(0);
    const int slice_extent_1 = precip_frac_2.extent_int(1);
    Kokkos::parallel_for(
        "precip_fraction::array_3", slice_extent_0 * slice_extent_1, KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = cell / (1) % slice_extent_0 + 1;
          const int slice_1 = cell / (1 * slice_extent_0) % slice_extent_1 + 1;
          precip_frac_2(slice_0 - 1, slice_1 - 1) = zero;
        });
  }
  // Set the minimum allowable precipitation fraction when hydrometeors are
  // found at a grid level.
  if (frozen_hm_present) {
    // Ice microphysics included.
    {
      const int extent_0 = ((ngrdcol) - (1)) / (1) + 1;
      Kokkos::parallel_for(
          "precip_fraction::cells_4", extent_0, KOKKOS_LAMBDA(const int cell) {
            const int j = (1) + (cell / (1) % extent_0) * (1);
            real_t cloud_frac_maxval = -std::numeric_limits<real_t>::max();
            for (int level = 1; level <= nzt; ++level) {
              cloud_frac_maxval = cloud_frac(j - 1, level - 1) > cloud_frac_maxval
                                      ? cloud_frac(j - 1, level - 1)
                                      : cloud_frac_maxval;
            }
            real_t ice_supersat_frac_maxval = -std::numeric_limits<real_t>::max();
            for (int level = 1; level <= nzt; ++level) {
              ice_supersat_frac_maxval =
                  ice_supersat_frac(j - 1, level - 1) > ice_supersat_frac_maxval
                      ? ice_supersat_frac(j - 1, level - 1)
                      : ice_supersat_frac_maxval;
            }
            precip_frac_tol(j - 1) =
                ((precip_frac_tol_coef * (cloud_frac_maxval > ice_supersat_frac_maxval
                                              ? cloud_frac_maxval
                                              : ice_supersat_frac_maxval)) > cloud_frac_min
                     ? (precip_frac_tol_coef * (cloud_frac_maxval > ice_supersat_frac_maxval
                                                    ? cloud_frac_maxval
                                                    : ice_supersat_frac_maxval))
                     : cloud_frac_min);
          });
    }
  } else {
    // Warm microphysics.
    {
      const int extent_0 = ((ngrdcol) - (1)) / (1) + 1;
      Kokkos::parallel_for(
          "precip_fraction::cells_5", extent_0, KOKKOS_LAMBDA(const int cell) {
            const int j = (1) + (cell / (1) % extent_0) * (1);
            real_t cloud_frac_maxval = -std::numeric_limits<real_t>::max();
            for (int level = 1; level <= nzt; ++level) {
              cloud_frac_maxval = cloud_frac(j - 1, level - 1) > cloud_frac_maxval
                                      ? cloud_frac(j - 1, level - 1)
                                      : cloud_frac_maxval;
            }
            precip_frac_tol(j - 1) = ((precip_frac_tol_coef * cloud_frac_maxval) > cloud_frac_min
                                          ? (precip_frac_tol_coef * cloud_frac_maxval)
                                          : cloud_frac_min);
          });
    }
  }
  // !! Find overall precipitation fraction.
  // The precipitation fraction is the greatest cloud fraction at or above a
  // vertical level.
  if (frozen_hm_present) {
    // Ice microphysics included.
    {
      const int slice_extent_0 = precip_frac.extent_int(0);
      Kokkos::parallel_for(
          "precip_fraction::array_6", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            precip_frac(slice_0 - 1, gr_k_ub_zt - 1) =
                (cloud_frac(slice_0 - 1, gr_k_ub_zt - 1) >
                         ice_supersat_frac(slice_0 - 1, gr_k_ub_zt - 1)
                     ? cloud_frac(slice_0 - 1, gr_k_ub_zt - 1)
                     : ice_supersat_frac(slice_0 - 1, gr_k_ub_zt - 1));
          });
    }
    // Independent columns retain Native physical top-to-bottom recurrence.
    Kokkos::parallel_for(
        "precip_fraction::recurrence_7", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int j = column + 1;
          for (int k = (gr_k_ub_zt - gr_grid_dir_indx);
               ((-(gr_grid_dir_indx)) > 0 ? k <= gr_k_lb_zt : k >= gr_k_lb_zt);
               k += -(gr_grid_dir_indx)) {
            precip_frac(j - 1, k - 1) =
                ((precip_frac(j - 1, (k + gr_grid_dir_indx) - 1) > cloud_frac(j - 1, k - 1)
                      ? precip_frac(j - 1, (k + gr_grid_dir_indx) - 1)
                      : cloud_frac(j - 1, k - 1)) > ice_supersat_frac(j - 1, k - 1)
                     ? (precip_frac(j - 1, (k + gr_grid_dir_indx) - 1) > cloud_frac(j - 1, k - 1)
                            ? precip_frac(j - 1, (k + gr_grid_dir_indx) - 1)
                            : cloud_frac(j - 1, k - 1))
                     : ice_supersat_frac(j - 1, k - 1));
          }
        });
  } else {
    // Warm microphysics.
    {
      const int slice_extent_0 = precip_frac.extent_int(0);
      Kokkos::parallel_for(
          "precip_fraction::array_8", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            precip_frac(slice_0 - 1, gr_k_ub_zt - 1) = cloud_frac(slice_0 - 1, gr_k_ub_zt - 1);
          });
    }
    // Independent columns retain Native physical top-to-bottom recurrence.
    Kokkos::parallel_for(
        "precip_fraction::recurrence_9", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int j = column + 1;
          for (int k = (gr_k_ub_zt - gr_grid_dir_indx);
               ((-(gr_grid_dir_indx)) > 0 ? k <= gr_k_lb_zt : k >= gr_k_lb_zt);
               k += -(gr_grid_dir_indx)) {
            precip_frac(j - 1, k - 1) =
                (precip_frac(j - 1, (k + gr_grid_dir_indx) - 1) > cloud_frac(j - 1, k - 1)
                     ? precip_frac(j - 1, (k + gr_grid_dir_indx) - 1)
                     : cloud_frac(j - 1, k - 1));
          }
        });
  }
  // !! Special checks for overall precipitation fraction
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "precip_fraction::cells_10", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          bool any_hydromet = false;
          bool all_hydromet_below = true;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(j - 1, k - 1, species - 1) >= hydromet_tol(species - 1));
            all_hydromet_below = all_hydromet_below &&
                                 (hydromet(j - 1, k - 1, species - 1) < hydromet_tol(species - 1));
          }
          if ((any_hydromet && (precip_frac(j - 1, k - 1) < precip_frac_tol(j - 1)))) {
            // In a scenario where we find any hydrometeor at this grid level, but
            // no cloud at or above this grid level, set precipitation fraction to
            // a minimum threshold value.
            precip_frac(j - 1, k - 1) = precip_frac_tol(j - 1);
          } else if (all_hydromet_below) {
            // The means (overall) of every precipitating hydrometeor are all less
            // than their respective tolerance amounts.  They are all considered to
            // have values of 0.  There are not any hydrometeor species found at
            // this grid level.  There is also no cloud at or above this grid
            // level, so set precipitation fraction to 0.
            precip_frac(j - 1, k - 1) = zero;
          }
        });
  }
  // Special checks for overall precipitation fraction loop: k = 1, nzt, 1
  // !! Find precipitation fraction within each PDF component.
  // The overall precipitation fraction, f_p, is given by the equation:
  // f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
  // where "a" is the mixture fraction (weight of PDF component 1), f_p(1) is
  // the precipitation fraction within PDF component 1, and f_p(2) is the
  // precipitation fraction within PDF component 2.  Overall precipitation
  // fraction is found according the method above, and mixture fraction is
  // already determined, leaving f_p(1) and f_p(2) to be solved for.  The
  // values for f_p(1) and f_p(2) must satisfy the above equation.
  if ((precip_frac_calc_type == 1)) {
    // Calculatate precip_frac_1 and precip_frac_2 based on the greatest
    // weighted cloud_frac_1 at or above a grid level.
    component_precip_frac_weighted(gr, nzt, ngrdcol, hydromet_dim, // In
                                   l_frozen_hm, hydromet_tol,      // In
                                   hydromet, precip_frac,          // In
                                   cloud_frac_1, cloud_frac_2,     // In
                                   ice_supersat_frac_1,            // In
                                   ice_supersat_frac_2, mixt_frac, // In
                                   precip_frac_tol,                // In
                                   precip_frac_1, precip_frac_2    // Out
    );
  } else if ((precip_frac_calc_type == 2)) {
    // Specified method.
    component_precip_frac_specify(
        nzt, ngrdcol, hydromet_dim, hydromet_tol,                             // In
        Kokkos::subview(clubb_params, Kokkos::ALL, iupsilon_precip_frac_rat), // In
        hydromet, precip_frac,                                                // In
        mixt_frac, precip_frac_tol,                                           // In
        precip_frac_1, precip_frac_2                                          // Out
    );
  } else {
    // Invalid option selected.
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "Invalid option to calculate precip_frac_1 and precip_frac_2." << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // precip_frac_calc_type
  // Increase Precipiation Fraction under special conditions.
  // There are scenarios that sometimes occur that require precipitation
  // fraction to be boosted.  Precipitation fraction is calculated from cloud
  // fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
  // PDF may become entirely subsaturated with respect to liquid and ice,
  // resulting in both a cloud fraction of 0 and an ice supersaturation
  // fraction of 0.  When this happens, precipitation fraction drops to 0 when
  // there aren't any hydrometeors present at that grid level, or to
  // precip_frac_tol when there is at least one hydrometeor present at that
  // grid level.  However, sometimes there are large values of hydrometeors
  // found at that grid level.  When this occurs, the PDF component in-precip
  // mean of a hydrometeor can become ridiculously large.  This is because the
  // ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
  // the equation:
  // mu_hm_i = hm_i / precip_frac_i;
  // where hm_i is the overall ith PDF component mean of the hydrometeor, and
  // precip_frac_i is the ith PDF component precipitation fraction.  When
  // precip_frac_i has a value of precip_frac_tol and hm_i is large, mu_hm_i
  // can be huge.  This can cause enormous microphysical process rates and
  // result in numerical instability.  It is also very inaccurate.
  // In order to limit this problem, the ith PDF component precipitation
  // fraction is increased in order to decrease mu_hm_i.  First, an "upper
  // limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
  // called max_hm_ip_comp_mean.  At every vertical level and for every
  // hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
  // exceeding the "upper limit".  The check is:
  // hm_i / precip_frac_i ( which = mu_hm_i )  >  max_hm_ip_comp_mean,
  // which can be rewritten:
  // hm_i > precip_frac_i * max_hm_ip_comp_mean.
  // Since hm_i has not been calculated yet, the assumption for this check is
  // that all of the hydrometeor is found in one PDF component, which is the
  // worst-case scenario in violating this limit.  The check becomes:
  // <hm> / ( mixt_frac * precip_frac_1 )  >  max_hm_ip_comp_mean;
  // in PDF comp. 1; and
  // <hm> / ( ( 1 - mixt_frac ) * precip_frac_2 )  >  max_hm_ip_comp_mean;
  // in PDF comp. 2.
  // These limits can be rewritten as:
  // <hm>  >  mixt_frac * precip_frac_1 * max_hm_ip_comp_mean;
  // in PDF comp. 1; and
  // <hm>  >  ( 1 - mixt_frac ) * precip_frac_2 * max_hm_ip_comp_mean;
  // in PDF comp. 2.
  // When component precipitation fraction is found to be in excess of the
  // limit, precip_frac_i is increased to:
  // <hm> / ( mixt_frac * max_hm_ip_comp_mean );
  // when the limit is exceeded in PDF comp. 1; and
  // <hm> / ( ( 1 - mixt_frac ) * max_hm_ip_comp_mean );
  // when the limit is exceeded in PDF comp. 2.
  // Of course, precip_frac_i is not allowed to exceed 1, so when
  // <hm> / mixt_frac (or <hm> / ( 1 - mixt_frac )) is already greater than
  // max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
  // max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
  // compared to what it would have been using precip_frac_tol.  In the event
  // that multiple hydrometeor mixing ratios violate the check, the code is set
  // up so that precip_frac_i is increased based on the highest hm_i.
  for (ivar = 1; ivar <= (hydromet_dim); ivar += 1) {
    if (l_mix_rat_hm_host(ivar - 1)) {
      {
        const int extent_0 = ((nzt) - (1)) / (1) + 1;
        const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
        Kokkos::parallel_for(
            "precip_fraction::cells_11", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int j = (1) + (cell / (1) % extent_1) * (1);
              const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
              // The hydrometeor is a mixing ratio.
              if (((hydromet(j - 1, k - 1, ivar - 1) >= hydromet_tol(ivar - 1)) &&
                   (hydromet(j - 1, k - 1, ivar - 1) >
                    ((mixt_frac(j - 1, k - 1) * precip_frac_1(j - 1, k - 1)) *
                     max_hm_ip_comp_mean)))) {
                // Increase precipitation fraction in the 1st PDF component.
                precip_frac_1(j - 1, k - 1) =
                    ((hydromet(j - 1, k - 1, ivar - 1) /
                      (mixt_frac(j - 1, k - 1) * max_hm_ip_comp_mean)) < one
                         ? (hydromet(j - 1, k - 1, ivar - 1) /
                            (mixt_frac(j - 1, k - 1) * max_hm_ip_comp_mean))
                         : one);
                // The value of precip_frac_1 must be at least precip_frac_tol
                // when precipitation is found in the 1st PDF component.
                precip_frac_1(j - 1, k - 1) = (precip_frac_1(j - 1, k - 1) > precip_frac_tol(j - 1)
                                                   ? precip_frac_1(j - 1, k - 1)
                                                   : precip_frac_tol(j - 1));
              }
              // <hm>/(mixt_frac*precip_frac_1) > max_hm_ip_comp_mean
              if (((hydromet(j - 1, k - 1, ivar - 1) >= hydromet_tol(ivar - 1)) &&
                   (hydromet(j - 1, k - 1, ivar - 1) >
                    (((one - mixt_frac(j - 1, k - 1)) * precip_frac_2(j - 1, k - 1)) *
                     max_hm_ip_comp_mean)))) {
                // Increase precipitation fraction in the 2nd PDF component.
                precip_frac_2(j - 1, k - 1) =
                    ((hydromet(j - 1, k - 1, ivar - 1) /
                      ((one - mixt_frac(j - 1, k - 1)) * max_hm_ip_comp_mean)) < one
                         ? (hydromet(j - 1, k - 1, ivar - 1) /
                            ((one - mixt_frac(j - 1, k - 1)) * max_hm_ip_comp_mean))
                         : one);
                // The value of precip_frac_2 must be at least precip_frac_tol
                // when precipitation is found in the 2nd PDF component.
                precip_frac_2(j - 1, k - 1) = (precip_frac_2(j - 1, k - 1) > precip_frac_tol(j - 1)
                                                   ? precip_frac_2(j - 1, k - 1)
                                                   : precip_frac_tol(j - 1));
              }
              // <hm>/((1-mixt_frac)*precip_frac_2) > max_hm_ip_comp_mean
            });
      }
      // k = 1, nzt
    }
    // l_mix_rat_hm(ivar)
  }
  // ivar = 1, hydromet_dim
  // Recalculate overall precipitation fraction for consistency.
  {
    const int slice_extent_0 = precip_frac.extent_int(0);
    const int slice_extent_1 = precip_frac.extent_int(1);
    Kokkos::parallel_for(
        "precip_fraction::array_12", slice_extent_0 * slice_extent_1,
        KOKKOS_LAMBDA(const int cell) {
          const int slice_0 = cell / (1) % slice_extent_0 + 1;
          const int slice_1 = cell / (1 * slice_extent_0) % slice_extent_1 + 1;
          precip_frac(slice_0 - 1, slice_1 - 1) =
              ((mixt_frac(slice_0 - 1, slice_1 - 1) * precip_frac_1(slice_0 - 1, slice_1 - 1)) +
               ((one - mixt_frac(slice_0 - 1, slice_1 - 1)) *
                precip_frac_2(slice_0 - 1, slice_1 - 1)));
        });
  }
  // Double check that precip_frac_tol <= precip_frac <= 1 when hydrometeors
  // are found at a grid level.
  // PLEASE DO NOT ALTER precip_frac, precip_frac_1, or precip_frac_2 anymore
  // after this point in the code.
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "precip_fraction::cells_13", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          bool any_hydromet = false;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(j - 1, k - 1, species - 1) >= hydromet_tol(species - 1));
          }
          if (any_hydromet) {
            precip_frac(j - 1, k - 1) =
                ((precip_frac(j - 1, k - 1) > precip_frac_tol(j - 1) ? precip_frac(j - 1, k - 1)
                                                                     : precip_frac_tol(j - 1)) < one
                     ? (precip_frac(j - 1, k - 1) > precip_frac_tol(j - 1)
                            ? precip_frac(j - 1, k - 1)
                            : precip_frac_tol(j - 1))
                     : one);
          }
          // any( hydromet(k,:) >= hydromet_tol(:) )
        });
  }
  // k = 1, nzt, 1
  // Statistics
  if (stats.flags().l_sample) {
    f2k::stats_update("precip_frac_tol", f2k::Real1D(precip_frac_tol), // In
                      stats                                            // InOut
    );
  }
  // Assertion check for precip_frac, precip_frac_1, and precip_frac_2.
  if (f2k::clubb_at_least_debug_level_api(2)) {
    const auto precip_frac_tol_host = to_host_space(precip_frac_tol);
    for (j = 1; j <= (ngrdcol); j += 1) {
      precip_frac_assert_check(nzt, hydromet_dim, hydromet_tol,                            // In
                               Kokkos::subview(hydromet, j - 1, Kokkos::ALL, Kokkos::ALL), // In
                               Kokkos::subview(mixt_frac, j - 1, Kokkos::ALL),             // In
                               Kokkos::subview(precip_frac, j - 1, Kokkos::ALL),           // In
                               Kokkos::subview(precip_frac_1, j - 1, Kokkos::ALL),         // In
                               Kokkos::subview(precip_frac_2, j - 1, Kokkos::ALL),         // In
                               precip_frac_tol_host(j - 1),                                // In
                               err_info                                                    // InOut
      );
      if ((to_host_space(err_info.err_code)(j - 1) == clubb_fatal_error)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          const auto header_host = to_host_space(err_info.err_header);
          diagnostic.write(header_host.data() + ((j)-1) * header_host.extent(0),
                           header_host.extent(0));
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "precip_frac_assert_check failed" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
    }
  }
  return;
}
// =============================================================================
namespace {
void component_precip_frac_weighted(const grid& gr, int nzt, int ngrdcol, int hydromet_dim, // In
                                    Bool1D l_frozen_hm, Real1D hydromet_tol,                // In
                                    Real3D hydromet, Real2D precip_frac,                    // In
                                    Real2D cloud_frac_1, Real2D cloud_frac_2,               // In
                                    Real2D ice_supersat_frac_1,                             // In
                                    Real2D ice_supersat_frac_2, Real2D mixt_frac,           // In
                                    Real1D precip_frac_tol,                                 // In
                                    Real2D precip_frac_1, Real2D precip_frac_2              // Out
) {
  // Description:
  // Set precipitation fraction in each component of the PDF.  The weighted 1st
  // PDF component precipitation fraction (weighted_pfrac_1) at a grid level is
  // calculated by the greatest value of mixt_frac * cloud_frac_1 at or above
  // the relevant grid level.  Likewise, the weighted 2nd PDF component
  // precipitation fraction (weighted_pfrac_2) at a grid level is calculated by
  // the greatest value of ( 1 - mixt_frac ) * cloud_frac_2 at or above the
  // relevant grid level.
  // The fraction weighted_pfrac_1 / ( weighted_pfrac_1 + weighted_pfrac_2 ) is
  // the weighted_pfrac_1 fraction.  Multiplying this fraction by overall
  // precipitation fraction and then dividing by mixt_frac produces the 1st PDF
  // component precipitation fraction (precip_frac_1).  Then, calculate the 2nd
  // PDF component precipitation fraction (precip_frac_2) accordingly.
  // References:
  //-----------------------------------------------------------------------

  // gr: Grid type variable
  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // l_frozen_hm: if true, then the hydrometeor is frozen; otherwise liquid
  // hydromet_tol: Tolerance values for all hydrometeors    [units vary]
  // hydromet: Mean of hydrometeor, hm (overall)           [units vary]
  // precip_frac: Precipitation fraction (overall)              [-]
  // cloud_frac_1: Cloud fraction (1st PDF component)            [-]
  // cloud_frac_2: Cloud fraction (2nd PDF component)            [-]
  // ice_supersat_frac_1: Ice supersaturation fraction (1st PDF comp.)  [-]
  // ice_supersat_frac_2: Ice supersaturation fraction (2nd PDF comp.)  [-]
  // mixt_frac: Mixture fraction                              [-]
  // precip_frac_tol: Minimum precip. frac. when hydromet. are present  [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)     [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)     [-]

  Real2D weighted_pfrac_1("component_precip_frac_weighted::weighted_pfrac_1", ngrdcol,
                          nzt); // Product of mixt_frac and cloud_frac_1           [-]
  Real2D weighted_pfrac_2("component_precip_frac_weighted::weighted_pfrac_2", ngrdcol,
                          nzt); // Product of ( 1 - mixt_frac ) and cloud_frac_2   [-]
  // Reduce immutable species flags once per call; no floating reduction.
  int frozen_hm_count = 0;
  Kokkos::parallel_reduce(
      "component_precip_frac_weighted::frozen_species", hydromet_dim,
      KOKKOS_LAMBDA(const int species, int& count) { count += l_frozen_hm(species) != 0; },
      frozen_hm_count);
  const bool frozen_hm_present = frozen_hm_count != 0;
  // Capture component Views/scalars only; owning UDTs remain on the host.
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zt = gr.k_lb_zt;

  // !! Find precipitation fraction within PDF component 1.
  // The method used to find overall precipitation fraction will also be to
  // find precipitation fraction within PDF component 1 and PDF component 2.
  // In order to do so, it is assumed (poorly) that PDF component 1 overlaps
  // PDF component 1 and that PDF component 2 overlaps PDF component 2 at every
  // vertical level in the vertical profile.
  // The weighted precipitation fraction in PDF component 1 is the greatest
  // value of the product of mixture fraction (mixt_frac) and 1st PDF
  // component cloud fraction at or above a vertical level.  Likewise, the
  // weighted precipitation fraction in PDF component 2 is the greatest
  // value of the product of ( 1 - mixt_frac ) and 2nd PDF component cloud
  // fraction at or above a vertical level.
  if (frozen_hm_present) {
    // Ice microphysics included.
    // Weighted precipitation fraction in PDF component 1.
    {
      const int slice_extent_0 = weighted_pfrac_1.extent_int(0);
      Kokkos::parallel_for(
          "component_precip_frac_weighted::array_1", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            weighted_pfrac_1(slice_0 - 1, gr_k_ub_zt - 1) =
                ((mixt_frac(slice_0 - 1, gr_k_ub_zt - 1) *
                  cloud_frac_1(slice_0 - 1, gr_k_ub_zt - 1)) >
                         (mixt_frac(slice_0 - 1, gr_k_ub_zt - 1) *
                          ice_supersat_frac_1(slice_0 - 1, gr_k_ub_zt - 1))
                     ? (mixt_frac(slice_0 - 1, gr_k_ub_zt - 1) *
                        cloud_frac_1(slice_0 - 1, gr_k_ub_zt - 1))
                     : (mixt_frac(slice_0 - 1, gr_k_ub_zt - 1) *
                        ice_supersat_frac_1(slice_0 - 1, gr_k_ub_zt - 1)));
          });
    }
    // Weighted precipitation fraction in PDF component 2.
    {
      const int slice_extent_0 = weighted_pfrac_2.extent_int(0);
      Kokkos::parallel_for(
          "component_precip_frac_weighted::array_2", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            weighted_pfrac_2(slice_0 - 1, gr_k_ub_zt - 1) =
                (((one - mixt_frac(slice_0 - 1, gr_k_ub_zt - 1)) *
                  cloud_frac_2(slice_0 - 1, gr_k_ub_zt - 1)) >
                         ((one - mixt_frac(slice_0 - 1, gr_k_ub_zt - 1)) *
                          ice_supersat_frac_2(slice_0 - 1, gr_k_ub_zt - 1))
                     ? ((one - mixt_frac(slice_0 - 1, gr_k_ub_zt - 1)) *
                        cloud_frac_2(slice_0 - 1, gr_k_ub_zt - 1))
                     : ((one - mixt_frac(slice_0 - 1, gr_k_ub_zt - 1)) *
                        ice_supersat_frac_2(slice_0 - 1, gr_k_ub_zt - 1)));
          });
    }
    // Independent columns retain Native physical top-to-bottom recurrence.
    Kokkos::parallel_for(
        "component_precip_frac_weighted::recurrence_3", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int j = column + 1;
          for (int k = (gr_k_ub_zt - gr_grid_dir_indx);
               ((-(gr_grid_dir_indx)) > 0 ? k <= gr_k_lb_zt : k >= gr_k_lb_zt);
               k += -(gr_grid_dir_indx)) {
            // Weighted precipitation fraction in PDF component 1.
            weighted_pfrac_1(j - 1, k - 1) =
                ((weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1) >
                          (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1))
                      ? weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1)
                      : (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1))) >
                         (mixt_frac(j - 1, k - 1) * ice_supersat_frac_1(j - 1, k - 1))
                     ? (weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1) >
                                (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1))
                            ? weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1)
                            : (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1)))
                     : (mixt_frac(j - 1, k - 1) * ice_supersat_frac_1(j - 1, k - 1)));
            // Weighted precipitation fraction in PDF component 2.
            weighted_pfrac_2(j - 1, k - 1) =
                ((weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1) >
                          ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1))
                      ? weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1)
                      : ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1))) >
                         ((one - mixt_frac(j - 1, k - 1)) * ice_supersat_frac_2(j - 1, k - 1))
                     ? (weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1) >
                                ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1))
                            ? weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1)
                            : ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1)))
                     : ((one - mixt_frac(j - 1, k - 1)) * ice_supersat_frac_2(j - 1, k - 1)));
          }
        });
  } else {
    // Warm microphysics.
    // Weighted precipitation fraction in PDF component 1.
    {
      const int slice_extent_0 = weighted_pfrac_1.extent_int(0);
      Kokkos::parallel_for(
          "component_precip_frac_weighted::array_4", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            weighted_pfrac_1(slice_0 - 1, gr_k_ub_zt - 1) =
                (mixt_frac(slice_0 - 1, gr_k_ub_zt - 1) *
                 cloud_frac_1(slice_0 - 1, gr_k_ub_zt - 1));
          });
    }
    // Weighted precipitation fraction in PDF component 2.
    {
      const int slice_extent_0 = weighted_pfrac_2.extent_int(0);
      Kokkos::parallel_for(
          "component_precip_frac_weighted::array_5", slice_extent_0, KOKKOS_LAMBDA(const int cell) {
            const int slice_0 = cell / (1) % slice_extent_0 + 1;
            weighted_pfrac_2(slice_0 - 1, gr_k_ub_zt - 1) =
                ((one - mixt_frac(slice_0 - 1, gr_k_ub_zt - 1)) *
                 cloud_frac_2(slice_0 - 1, gr_k_ub_zt - 1));
          });
    }
    // Independent columns retain Native physical top-to-bottom recurrence.
    Kokkos::parallel_for(
        "component_precip_frac_weighted::recurrence_6", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int j = column + 1;
          for (int k = (gr_k_ub_zt - gr_grid_dir_indx);
               ((-(gr_grid_dir_indx)) > 0 ? k <= gr_k_lb_zt : k >= gr_k_lb_zt);
               k += -(gr_grid_dir_indx)) {
            // Weighted precipitation fraction in PDF component 1.
            weighted_pfrac_1(j - 1, k - 1) =
                (weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1) >
                         (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1))
                     ? weighted_pfrac_1(j - 1, (k + gr_grid_dir_indx) - 1)
                     : (mixt_frac(j - 1, k - 1) * cloud_frac_1(j - 1, k - 1)));
            // Weighted precipitation fraction in PDF component 2.
            weighted_pfrac_2(j - 1, k - 1) =
                (weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1) >
                         ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1))
                     ? weighted_pfrac_2(j - 1, (k + gr_grid_dir_indx) - 1)
                     : ((one - mixt_frac(j - 1, k - 1)) * cloud_frac_2(j - 1, k - 1)));
          }
        });
  }
  // Calculate precip_frac_1 and special cases for precip_frac_1.
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "component_precip_frac_weighted::cells_7", extent_0 * extent_1,
        KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          // Calculate precipitation fraction in the 1st PDF component.
          if (((weighted_pfrac_1(j - 1, k - 1) + weighted_pfrac_2(j - 1, k - 1)) > zero)) {
            // Adjust weighted 1st PDF component precipitation fraction by
            // multiplying it by a factor.  That factor is overall precipitation
            // fraction divided by the sum of the weighted 1st PDF component
            // precipitation fraction and the weighted 2nd PDF component
            // precipitation fraction.  The 1st PDF component precipitation
            // fraction is then found by dividing the adjusted weighted 1st PDF
            // component precipitation fraction by mixture fraction.
            precip_frac_1(j - 1, k - 1) =
                ((weighted_pfrac_1(j - 1, k - 1) *
                  (precip_frac(j - 1, k - 1) /
                   (weighted_pfrac_1(j - 1, k - 1) + weighted_pfrac_2(j - 1, k - 1)))) /
                 mixt_frac(j - 1, k - 1));
          } else {
            // Usually, the sum of the weighted 1st PDF component precipitation
            // fraction and the 2nd PDF component precipitation fraction go to 0
            // when overall precipitation fraction goes to 0.  Since 1st PDF
            // component weighted precipitation fraction is 0, 1st PDF component
            // precipitation fraction also 0.
            precip_frac_1(j - 1, k - 1) = zero;
          }
        });
  }
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "component_precip_frac_weighted::cells_8", extent_0 * extent_1,
        KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          bool any_hydromet = false;
          bool all_hydromet_below = true;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(j - 1, k - 1, species - 1) >= hydromet_tol(species - 1));
            all_hydromet_below = all_hydromet_below &&
                                 (hydromet(j - 1, k - 1, species - 1) < hydromet_tol(species - 1));
          }
          // Special cases for precip_frac_1.
          if ((any_hydromet && (precip_frac_1(j - 1, k - 1) >
                                (one < (precip_frac(j - 1, k - 1) / mixt_frac(j - 1, k - 1))
                                     ? one
                                     : (precip_frac(j - 1, k - 1) / mixt_frac(j - 1, k - 1)))))) {
            // Using the above method, it is possible for precip_frac_1 to be
            // greater than 1.  For example, the mixture fraction at level k+1 is
            // 0.10 and the cloud_frac_1 at level k+1 is 1, resulting in a
            // weighted_pfrac_1 of 0.10.  This product is greater than the product
            // of mixt_frac and cloud_frac_1 at level k.  The mixture fraction at
            // level k is 0.05, resulting in a precip_frac_1 of 2.  The value of
            // precip_frac_1 is limited at 1.  The leftover precipitation fraction
            // (a result of the decreasing weight of PDF component 1 between the
            // levels) is applied to PDF component 2.
            // Additionally, when weighted_pfrac_1 at level k is greater than
            // overall precipitation fraction at level k, the resulting calculation
            // of precip_frac_2 at level k will be negative.
            precip_frac_1(j - 1, k - 1) =
                (one < (precip_frac(j - 1, k - 1) / mixt_frac(j - 1, k - 1))
                     ? one
                     : (precip_frac(j - 1, k - 1) / mixt_frac(j - 1, k - 1)));
          } else if ((any_hydromet && (precip_frac_1(j - 1, k - 1) > zero) &&
                      (precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1)))) {
            // In a scenario where we find precipitation in the 1st PDF component
            // (it is allowed to have a value of 0 when all precipitation is found
            // in the 2nd PDF component) but it is tiny (less than tolerance
            // level), boost 1st PDF component precipitation fraction to tolerance
            // level.
            precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
          } else if (all_hydromet_below) {
            // The means (overall) of every precipitating hydrometeor are all less
            // than their respective tolerance amounts.  They are all considered to
            // have values of 0.  There are not any hydrometeor species found at
            // this grid level.  There is also no cloud at or above this grid
            // level, so set 1st component precipitation fraction to 0.
            precip_frac_1(j - 1, k - 1) = zero;
          }
        });
  }
  // !! Find precipitation fraction within PDF component 2.
  // The equation for precipitation fraction within PDF component 2 is:
  // f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
  // given the overall precipitation fraction, f_p (calculated above), the
  // precipitation fraction within PDF component 1, f_p(1) (calculated above),
  // and mixture fraction, a.  Any leftover precipitation fraction from
  // precip_frac_1 will be included in this calculation of precip_frac_2.
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "component_precip_frac_weighted::cells_9", extent_0 * extent_1,
        KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          bool any_hydromet = false;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(j - 1, k - 1, species - 1) >= hydromet_tol(species - 1));
          }
          if (any_hydromet) {
            // Calculate precipitation fraction in the 2nd PDF component.
            precip_frac_2(j - 1, k - 1) =
                (((precip_frac(j - 1, k - 1) -
                   (mixt_frac(j - 1, k - 1) * precip_frac_1(j - 1, k - 1))) /
                  (one - mixt_frac(j - 1, k - 1))) > zero
                     ? ((precip_frac(j - 1, k - 1) -
                         (mixt_frac(j - 1, k - 1) * precip_frac_1(j - 1, k - 1))) /
                        (one - mixt_frac(j - 1, k - 1)))
                     : zero);
            // Special cases for precip_frac_2.
            if ((precip_frac_2(j - 1, k - 1) > one)) {
              // Again, it is possible for precip_frac_2 to be greater than 1.
              // For example, the mixture fraction at level k+1 is 0.10 and the
              // cloud_frac_1 at level k+1 is 1, resulting in a weighted_pfrac_1
              // of 0.10.  This product is greater than the product of mixt_frac
              // and cloud_frac_1 at level k.  Additionally, precip_frac (overall)
              // is 1 for level k.  The mixture fraction at level k is 0.5,
              // resulting in a precip_frac_1 of 0.2.  Using the above equation,
              // precip_frac_2 is calculated to be 1.8.  The value of
              // precip_frac_2 is limited at 1.  The leftover precipitation
              // fraction (as a result of the increasing weight of component 1
              // between the levels) is applied to PDF component 1.
              precip_frac_2(j - 1, k - 1) = one;
              // Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j - 1, k - 1) =
                  ((precip_frac(j - 1, k - 1) - (one - mixt_frac(j - 1, k - 1))) /
                   mixt_frac(j - 1, k - 1));
              // Double check precip_frac_1
              if ((precip_frac_1(j - 1, k - 1) > one)) {
                precip_frac_1(j - 1, k - 1) = one;
                precip_frac_2(j - 1, k - 1) =
                    ((precip_frac(j - 1, k - 1) - mixt_frac(j - 1, k - 1)) /
                     (one - mixt_frac(j - 1, k - 1)));
              } else if (((precip_frac_1(j - 1, k - 1) > zero) &&
                          (precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1)))) {
                precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
                // fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j - 1, k - 1) =
                    (precip_frac_1(j - 1, k - 1) *
                     (((precip_frac(j - 1, k - 1) / precip_frac_1(j - 1, k - 1)) -
                       mixt_frac(j - 1, k - 1)) /
                      (one - mixt_frac(j - 1, k - 1))));
              }
            } else if (((precip_frac_2(j - 1, k - 1) > zero) &&
                        (precip_frac_2(j - 1, k - 1) < precip_frac_tol(j - 1)))) {
              // In a scenario where we find precipitation in the 2nd PDF
              // component (it is allowed to have a value of 0 when all
              // precipitation is found in the 1st PDF component) but it is tiny
              // (less than tolerance level), boost 2nd PDF component
              // precipitation fraction to tolerance level.
              precip_frac_2(j - 1, k - 1) = precip_frac_tol(j - 1);
              // Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j - 1, k - 1) =
                  ((precip_frac(j - 1, k - 1) -
                    ((one - mixt_frac(j - 1, k - 1)) * precip_frac_2(j - 1, k - 1))) /
                   mixt_frac(j - 1, k - 1));
              // Double check precip_frac_1
              if ((precip_frac_1(j - 1, k - 1) > one)) {
                precip_frac_1(j - 1, k - 1) = one;
                precip_frac_2(j - 1, k - 1) =
                    ((precip_frac(j - 1, k - 1) - mixt_frac(j - 1, k - 1)) /
                     (one - mixt_frac(j - 1, k - 1)));
              } else if (((precip_frac_1(j - 1, k - 1) > zero) &&
                          (precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1)))) {
                precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
                // fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j - 1, k - 1) =
                    (precip_frac_1(j - 1, k - 1) *
                     (((precip_frac(j - 1, k - 1) / precip_frac_1(j - 1, k - 1)) -
                       mixt_frac(j - 1, k - 1)) /
                      (one - mixt_frac(j - 1, k - 1))));
              }
            }
            // Special cases for precip_frac_2
          } else {
            // all( hydromet(k,:) < hydromet_tol(:) )
            // The means (overall) of every precipitating hydrometeor are all less
            // than their respective tolerance amounts.  They are all considered to
            // have values of 0.  There are not any hydrometeor species found at
            // this grid level.  There is also no cloud at or above this grid
            // level, so set 2nd component precipitation fraction to 0.
            precip_frac_2(j - 1, k - 1) = zero;
          }
          // any( hydromet(k,:) > hydromet_tol(:) )
        });
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void component_precip_frac_specify(int nzt, int ngrdcol, int hydromet_dim,
                                   Real1D hydromet_tol,                       // In
                                   Real1D upsilon_precip_frac_rat,            // In
                                   Real3D hydromet, Real2D precip_frac,       // In
                                   Real2D mixt_frac, Real1D precip_frac_tol,  // In
                                   Real2D precip_frac_1, Real2D precip_frac_2 // Out
) {
  // Description:
  // Calculates the precipitation fraction in each PDF component.
  // The equation for precipitation fraction is:
  // f_p = mixt_frac * f_p(1) + ( 1 - mixt_frac ) * f_p(2);
  // where f_p is overall precipitation fraction, f_p(1) is precipitation
  // fraction in the 1st PDF component, f_p(2) is precipitation fraction in the
  // 2nd PDF component, and mixt_frac is the mixture fraction.  Using this
  // method, a new specified parameter is introduced, upsilon, where:
  // upsilon = mixt_frac * f_p(1) / f_p; and where 0 <= upsilon <= 1.
  // In other words, upsilon is the ratio of mixt_frac * f_p(1) to f_p.  Since
  // f_p and mixt_frac are calculated previously, and upsilon is specified,
  // f_p(1) can be calculated by:
  // f_p(1) = upsilon * f_p / mixt_frac;
  // and has an upper limit of 1.  The value of f_p(2) can then be calculated
  // by:
  // f_p(2) = ( f_p - mixt_frac * f_p(1) ) / ( 1 - mixt_frac );
  // and also has an upper limit of 1.  When upsilon = 1, almost all of the
  // precipitation is found in the 1st PDF component (as long as
  // f_p <= mixt_frac, otherwise it would cause f_p(1) to be greater than 1).
  // When upsilon = 0, almost all of the precipitation is found in the 2nd PDF
  // component (as long as f_p <= 1 - mixt_frac, otherwise it would cause
  // f_p(2) to be greater than 1). Precipitation is split between the two PDF
  // components according to upsilon.
  // References:
  //-----------------------------------------------------------------------

  // nzt: Number of model thermodynamic vertical grid levels
  // ngrdcol: Number of grid columns
  // hydromet_dim: Number of hydrometeor species
  // hydromet_tol: Tolerance values for all hydrometeors    [units vary]
  // upsilon_precip_frac_rat: ratio mixt_frac*precip_frac_1/precip_frac [-]
  // hydromet: Mean of hydrometeor, hm (overall)           [units vary]
  // precip_frac: Precipitation fraction (overall)                      [-]
  // mixt_frac: Mixture fraction                                      [-]
  // precip_frac_tol: Minimum precip. frac. when hydromet. are present  [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)     [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)     [-]


  // Loop over all vertical levels.
  {
    const int extent_0 = ((nzt) - (1)) / (1) + 1;
    const int extent_1 = ((ngrdcol) - (1)) / (1) + 1;
    Kokkos::parallel_for(
        "component_precip_frac_specify::cells_1", extent_0 * extent_1,
        KOKKOS_LAMBDA(const int cell) {
          const int j = (1) + (cell / (1) % extent_1) * (1);
          const int k = (1) + (cell / (1 * extent_1) % extent_0) * (1);
          bool any_hydromet = false;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(j - 1, k - 1, species - 1) >= hydromet_tol(species - 1));
          }
          if (any_hydromet) {
            // Precipitation is found in both PDF components.  Each component
            // must have a precipitation fraction that is at least
            // precip_frac_tol and that does not exceed 1.
            // Calculate precipitation fraction in the 1st PDF component.
            precip_frac_1(j - 1, k - 1) =
                ((upsilon_precip_frac_rat(j - 1) * precip_frac(j - 1, k - 1)) /
                 mixt_frac(j - 1, k - 1));
            // Special cases for precip_frac_1
            if ((precip_frac_1(j - 1, k - 1) > one)) {
              precip_frac_1(j - 1, k - 1) = one;
            } else if ((precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1))) {
              precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
            }
            // Calculate precipitation fraction in the 2nd PDF component.
            precip_frac_2(j - 1, k - 1) =
                ((precip_frac(j - 1, k - 1) -
                  (mixt_frac(j - 1, k - 1) * precip_frac_1(j - 1, k - 1))) /
                 (one - mixt_frac(j - 1, k - 1)));
            // Special case for precip_frac_2
            if ((precip_frac_2(j - 1, k - 1) > one)) {
              // Set precip_frac_2 to 1.
              precip_frac_2(j - 1, k - 1) = one;
              // Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j - 1, k - 1) =
                  ((precip_frac(j - 1, k - 1) - (one - mixt_frac(j - 1, k - 1))) /
                   mixt_frac(j - 1, k - 1));
              // Double check precip_frac_1
              if ((precip_frac_1(j - 1, k - 1) > one)) {
                precip_frac_1(j - 1, k - 1) = one;
                precip_frac_2(j - 1, k - 1) =
                    ((precip_frac(j - 1, k - 1) - mixt_frac(j - 1, k - 1)) /
                     (one - mixt_frac(j - 1, k - 1)));
              } else if ((precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1))) {
                precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
                // fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j - 1, k - 1) =
                    (precip_frac_1(j - 1, k - 1) *
                     (((precip_frac(j - 1, k - 1) / precip_frac_1(j - 1, k - 1)) -
                       mixt_frac(j - 1, k - 1)) /
                      (one - mixt_frac(j - 1, k - 1))));
              }
            } else if ((precip_frac_2(j - 1, k - 1) < precip_frac_tol(j - 1))) {
              // Set precip_frac_2 to precip_frac_tol.
              precip_frac_2(j - 1, k - 1) = precip_frac_tol(j - 1);
              // Recalculate precipitation fraction in the 1st PDF component.
              precip_frac_1(j - 1, k - 1) =
                  ((precip_frac(j - 1, k - 1) -
                    ((one - mixt_frac(j - 1, k - 1)) * precip_frac_2(j - 1, k - 1))) /
                   mixt_frac(j - 1, k - 1));
              // Double check precip_frac_1
              if ((precip_frac_1(j - 1, k - 1) > one)) {
                precip_frac_1(j - 1, k - 1) = one;
                precip_frac_2(j - 1, k - 1) =
                    ((precip_frac(j - 1, k - 1) - mixt_frac(j - 1, k - 1)) /
                     (one - mixt_frac(j - 1, k - 1)));
              } else if ((precip_frac_1(j - 1, k - 1) < precip_frac_tol(j - 1))) {
                precip_frac_1(j - 1, k - 1) = precip_frac_tol(j - 1);
                // fp = a*fp1+(1-a)*fp2 solving for fp2
                precip_frac_2(j - 1, k - 1) =
                    (precip_frac_1(j - 1, k - 1) *
                     (((precip_frac(j - 1, k - 1) / precip_frac_1(j - 1, k - 1)) -
                       mixt_frac(j - 1, k - 1)) /
                      (one - mixt_frac(j - 1, k - 1))));
              }
            }
            // Special cases for precip_frac_2
          } else {
            // all( hydromet(k,:) < hydromet_tol(:) )
            // There aren't any hydrometeors found at the grid level.
            precip_frac_1(j - 1, k - 1) = zero;
            precip_frac_2(j - 1, k - 1) = zero;
          }
        });
  }
  return;
}
} // namespace
// =============================================================================
namespace {
void precip_frac_assert_check(int nzt, int hydromet_dim, Real1D hydromet_tol,        // In
                              Real2D hydromet, Real1D mixt_frac, Real1D precip_frac, // In
                              Real1D precip_frac_1, Real1D precip_frac_2,            // In
                              real_t precip_frac_tol,                                // In
                              err_info_type& err_info                                // InOut
) {
  // Description: Assertion check for precipitation fractions. Native first-failed
  // level/check order, one-sided consistency predicate and all-column fatal
  // assignment are retained. The device returns only failure metadata/values.
  // nzt: Number of model thermodynamic vertical grid levels
  // hydromet_dim: Number of hydrometeor species
  // hydromet_tol: Tolerance values for all hydrometeors    [units vary]
  // hydromet: Mean of hydrometeor, hm (overall)           [units vary]
  // mixt_frac: Mixture fraction                               [-]
  // precip_frac: Precipitation fraction (overall)               [-]
  // precip_frac_1: Precipitation fraction (1st PDF component)     [-]
  // precip_frac_2: Precipitation fraction (2nd PDF component)     [-]
  // precip_frac_tol: Minimum precip. frac. when hydromet. are present  [-]
  // err_info: err_info struct containing err_code and err_header

  Kokkos::View<int*, Kokkos::LayoutLeft> failure_code("precip_frac_assert_check::failure", 2);
  Real1D failure_values("precip_frac_assert_check::values", 4);
  Kokkos::parallel_for(
      "precip_frac_assert_check", 1, KOKKOS_LAMBDA(const int) {
        failure_code(0) = 0;
        failure_code(1) = 0;
        // Loop over all vertical levels.
        for (int k = 1; k <= nzt; k += 1) {
          bool any_hydromet = false;
          for (int species = 1; species <= hydromet_dim; ++species) {
            any_hydromet =
                any_hydromet || (hydromet(k - 1, species - 1) >= hydromet_tol(species - 1));
          }
          if (any_hydromet) {
            // Overall precipitation fraction cannot be less than precip_frac_tol
            // when a hydrometeor is present at a grid level.
            if ((precip_frac(k - 1) < precip_frac_tol)) {
              failure_code(0) = 1;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Overall precipitation fraction cannot exceed 1.
            if ((precip_frac(k - 1) > one)) {
              failure_code(0) = 2;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 1st PDF component is allowed to be 0
            // when all the precipitation is found in the 2nd PDF component.
            // Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
            // is present at a grid level.  In other words, it cannot have a value
            // that is greater than 0 but less than precip_frac_tol
            if (((precip_frac_1(k - 1) > zero) &&
                 (precip_frac_1(k - 1) < (precip_frac_tol - eps)))) {
              failure_code(0) = 3;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 1st PDF component cannot exceed 1.
            if ((precip_frac_1(k - 1) > one)) {
              failure_code(0) = 4;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipiation fraction in the 1st PDF component cannot be negative.
            if ((precip_frac_1(k - 1) < zero)) {
              failure_code(0) = 5;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 2nd PDF component is allowed to be 0
            // when all the precipitation is found in the 1st PDF component.
            // Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
            // is present at a grid level.  In other words, it cannot have a value
            // that is greater than 0 but less than precip_frac_tol
            if (((precip_frac_2(k - 1) > zero) &&
                 (precip_frac_2(k - 1) < (precip_frac_tol - eps)))) {
              failure_code(0) = 6;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 2nd PDF component cannot exceed 1.
            if ((precip_frac_2(k - 1) > one)) {
              failure_code(0) = 7;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipiation fraction in the 2nd PDF component cannot be negative.
            if ((precip_frac_2(k - 1) < zero)) {
              failure_code(0) = 8;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
          } else {
            // all( hydromet(k,:) < hydromet_tol(:) )
            // Overall precipitation fraction must be 0 when no hydrometeors are
            // found at a grid level.
            if ((Kokkos::abs(precip_frac(k - 1)) > eps)) {
              failure_code(0) = 9;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 1st PDF component must be 0 when no
            // hydrometeors are found at a grid level.
            if ((Kokkos::abs(precip_frac_1(k - 1)) > eps)) {
              failure_code(0) = 10;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
            // Precipitation fraction in the 2nd PDF component must be 0 when no
            // hydrometeors are found at a grid level.
            if ((Kokkos::abs(precip_frac_2(k - 1)) > eps)) {
              failure_code(0) = 11;
              failure_code(1) = k;
              failure_values(0) = mixt_frac(k - 1);
              failure_values(1) = precip_frac(k - 1);
              failure_values(2) = precip_frac_1(k - 1);
              failure_values(3) = precip_frac_2(k - 1);
              return;
            }
          }
          // any( hydromet(k,:) >= hydromet_tol(:) )
          // The precipitation fraction equation is:
          // precip_frac
          // = mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2;
          // which means that:
          // precip_frac
          // - ( mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2 )
          // = 0.
          // Check that this is true with numerical round off.
          if (((precip_frac(k - 1) - ((mixt_frac(k - 1) * precip_frac_1(k - 1)) +
                                      ((one - mixt_frac(k - 1)) * precip_frac_2(k - 1)))) > eps)) {
            failure_code(0) = 12;
            failure_code(1) = k;
            failure_values(0) = mixt_frac(k - 1);
            failure_values(1) = precip_frac(k - 1);
            failure_values(2) = precip_frac_1(k - 1);
            failure_values(3) = precip_frac_2(k - 1);
            return;
          }
        }
        // k = 1, nzt, 1
        return;
      });
  const auto failure_code_host = to_host_space(failure_code);
  if (failure_code_host(0) == 0)
    return;
  const int k = failure_code_host(1);
  const auto failure_values_host = to_host_space(failure_values);
  switch (failure_code_host(0)) {
  case 1: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac < precip_frac_tol when a hydrometeor is present" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac = " << ' ';
      diagnostic << failure_values_host(1) << ' ';
      diagnostic << "precip_frac_tol = " << ' ';
      diagnostic << precip_frac_tol << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 2: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac > 1" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac = " << ' ';
      diagnostic << failure_values_host(1) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 3: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "0 < precip_frac_1 < precip_frac_tol" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 = " << ' ';
      diagnostic << failure_values_host(2) << ' ';
      diagnostic << "precip_frac_tol = " << ' ';
      diagnostic << precip_frac_tol << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 4: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 > 1" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 = " << ' ';
      diagnostic << failure_values_host(2) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 5: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 < 0" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 = " << ' ';
      diagnostic << failure_values_host(2) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 6: {
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "0 < precip_frac_2 < precip_frac_tol" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 = " << ' ';
      diagnostic << failure_values_host(3) << ' ';
      diagnostic << "precip_frac_tol = " << ' ';
      diagnostic << precip_frac_tol << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 7: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 > 1" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 = " << ' ';
      diagnostic << failure_values_host(3) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 8: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 < 0" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 = " << ' ';
      diagnostic << failure_values_host(3) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 9: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac != 0 when no hydrometeors are found" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac = " << ' ';
      diagnostic << failure_values_host(1) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 10: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 /= 0 when no hydrometeors are found" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_1 = " << ' ';
      diagnostic << failure_values_host(2) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 11: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 /= 0 when no hydrometeors are found" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac_2 = " << ' ';
      diagnostic << failure_values_host(3) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  case 12: {
    std::cerr << err_info.err_header_global << "\n";
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2 /= precip_frac "
                    "within numerical roundoff"
                 << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "level = " << ' ';
      diagnostic << k << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2 = " << ' ';
      diagnostic << ((failure_values_host(0) * failure_values_host(2)) +
                     ((one - failure_values_host(0)) * failure_values_host(3)))
                 << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "precip_frac = " << ' ';
      diagnostic << failure_values_host(1) << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    break;
  }
  }
  // Native checker marks every column fatal, including when called for one
  // profile. Preserve that behavior rather than the old intended-column comment.
  Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
}
} // namespace
} // namespace clubb_kokkos
