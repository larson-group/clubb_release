// Advance-equation helpers on CLUBB's staggered thermodynamic/momentum grids.
//
// Derivatives and stability: calc_ddzt_umvm_sqd, calc_wp3_on_wp2,
// calc_stability_correction, calc_brunt_vaisala_freq_sqd, compute_Cx_fnc_Richardson
// and calc_Ri_zm diagnose shear, velocity ratios and stability corrections.
// Boundary conditions: set_boundary_conditions_lhs/rhs update LAPACK systems.
// Smoothing: Lscale_width_vert_avg supplies a density-weighted running mean;
// smooth_min/max select Native scalar/array formulas; smooth_heaviside_peskin
// smooths a step. wp23_term_splat_lhs diagnoses velocity-moment splatting.
// Other operations: calc_xpwp computes diffusive fluxes, vertical_avg and
// vertical_integral reduce profiles, pvertinterp interpolates pressure levels,
// and calculate_thlp2_rad updates radiative forcing of temperature variance.
//
// Arrays use column-first LayoutLeft execution-space Views. Generic boundary,
// smoothing and flux adapters also accept host/strided input; shared array_types
// conversions own transfers. Boundary-index arguments and grid metadata retain
// Native one-based values. Vertical reductions return host scalars and retain
// increasing-profile accumulation. Scratch is per call; stats stays Native.
// Scalar smooth_min/max formulas are device callable. No implicit bridge calls
// occur in kernels. See INTEGRATION.md for unintegrated draft status/dependencies.
// MIN/MAX ternaries retain the Native second operand on ties, as established
// by the grid port; integer powers use ordinary products. No explicit FMA or
// special-function workaround is introduced. Numerical evidence is pending.
// References and detailed scientific comments remain beside Native procedures.
#pragma once

#include "array_types.hpp"
#include "grid_class.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "T_in_K_module.hpp"
#include "saturation.hpp"
#include "f2k/stats_bridge.hpp"
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace clubb_kokkos {

// Forward declarations preserve Native definition order for lower calls.
inline void calc_ddzt_umvm_sqd(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D um,
                               Real2D vm,           // In
                               Real2D ddzt_umvm_sqd // Out
);
inline void calc_wp3_on_wp2(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D wp2,
                            Real2D wp3,                             // In
                            Real2D wp3_on_wp2, Real2D wp3_on_wp2_zt // Out
);
inline void
set_boundary_conditions_lhs(int diag_index, int low_bound, int high_bound, // In
                            Real2D lhs,                                    // InOut
                            std::optional<int> diag_index2 = std::nullopt, // Optional In
                            std::optional<int> low_bound2 = std::nullopt,  // Optional In
                            std::optional<int> high_bound2 = std::nullopt  // Optional In
);
inline void
set_boundary_conditions_rhs(real_t low_value, int low_bound, real_t high_value,
                            int high_bound,                                   // In
                            Real1D rhs,                                       // InOut
                            std::optional<real_t> low_value2 = std::nullopt,  // Optional In
                            std::optional<int> low_bound2 = std::nullopt,     // Optional In
                            std::optional<real_t> high_value2 = std::nullopt, // Optional In
                            std::optional<int> high_bound2 = std::nullopt     // Optional In
);
inline void calc_stability_correction(int nzm, int ngrdcol, // In
                                      Real2D brunt_vaisala_freq_sqd, Real2D Lscale_zm,
                                      Real2D em,                     // In
                                      Real1D lambda0_stability_coef, // In
                                      Real2D stability_correction    // Out
);
inline void calc_brunt_vaisala_freq_sqd(int nzm, int nzt, int ngrdcol, const grid& gr,
                                        Real2D thlm, // In
                                        Real2D exner, Real2D rtm, Real2D rcm, Real2D p_in_Pa,
                                        Real2D thvm,                          // In
                                        Real2D ice_supersat_frac,             // In
                                        int saturation_formula,               // In
                                        bool l_brunt_vaisala_freq_moist,      // In
                                        bool l_use_thvm_in_bv_freq,           // In
                                        bool l_modify_limiters_for_cnvg_test, // In
                                        Real1D bv_efold, real_t T0,           // In
                                        Real2D brunt_vaisala_freq_sqd,        // Out
                                        Real2D brunt_vaisala_freq_sqd_mixed,  // Out
                                        Real2D brunt_vaisala_freq_sqd_smth,   // Out
                                        f2k::stats_type* stats = nullptr      // Optional InOut
);
inline void compute_Cx_fnc_Richardson(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                      Real2D Lscale_zm, Real2D ddzt_umvm_sqd,
                                      Real2D rho_ds_zm,                     // In
                                      Real2D brunt_vaisala_freq_sqd,        // In
                                      Real2D brunt_vaisala_freq_sqd_mixed,  // In
                                      Real2D clubb_params,                  // In
                                      bool l_use_shear_Richardson,          // In
                                      bool l_modify_limiters_for_cnvg_test, // In
                                      Real2D Cx_fnc_Richardson              // Out
);
inline void calc_Ri_zm(int nzm, int ngrdcol, Real2D bv_freq_sqd, Real2D shear, real_t lim_bv, // In
                       real_t lim_shear,                                                      // In
                       Real2D Ri_zm                                                           // Out
);
inline Real2D Lscale_width_vert_avg(int nzm, int ngrdcol, const grid& gr, int smth_type,    // In
                                    Real2D var_profile, Real2D Lscale_zm, Real2D rho_ds_zm, // In
                                    real_t var_below_ground_value                           // In
);
inline void wp23_term_splat_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                                Real1D C_wp2_splat, // In
                                Real2D brunt_vaisala_freq_sqd_mixed, Real2D Lscale_zm,
                                Real2D rho_ds_zm,                          // In
                                f2k::stats_type& stats,                    // InOut
                                Real2D lhs_splat_wp2, Real2D lhs_splat_wp3 // Out
);
inline Real2D smooth_min_sclr_idx(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                                  real_t smth_coef // In
);
inline Real2D smooth_min_array_scalar(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                                      real_t smth_coef // In
);
inline Real2D smooth_min_arrays(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                                real_t smth_coef // In
);
KOKKOS_INLINE_FUNCTION real_t smooth_min_scalars(real_t input_var1, real_t input_var2,
                                                 real_t smth_coef // In
);
inline Real2D smooth_max_sclr_idx(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                                  real_t smth_coef // In
);
inline Real2D smooth_max_array_scalar(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                                      real_t smth_coef // In
);
inline Real1D smooth_max_array_1D_scalar(int ngrdcol, Real1D input_var1, real_t input_var2,
                                         real_t smth_coef // In
);
inline Real2D smooth_max_arrays(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                                real_t smth_coef // In
);
KOKKOS_INLINE_FUNCTION real_t smooth_max_scalars(real_t input_var1, real_t input_var2,
                                                 real_t smth_coef // In
);
inline Real2D smooth_heaviside_peskin(int nz, int ngrdcol, Real2D input, // In
                                      real_t smth_range                  // In
);
inline void calc_xpwp_1D(const grid& gr, Real1D Km_zm, Real1D xm, // In
                         Real1D xpwp                              // Out
);
inline void calc_xpwp_2D(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                         Real2D Km_zm, Real2D xm,                       // In
                         Real2D xpwp                                    // Out
);
inline void pvertinterp(int nzt, int ngrdcol, const grid& gr,         // In
                        Real2D p_mid, real_t p_out, Real2D input_var, // In
                        Real1D interp_var                             // Out
);
inline void calculate_thlp2_rad(int ngrdcol, int nzm, int nzt, const grid& gr,                 // In
                                Real2D rcm, Real2D thlprcp, Real2D radht, Real2D clubb_params, // In
                                Real2D thlp2_forcing // InOut
);

// C++ overload declarations implement the three Native generic interfaces.
inline Real2D smooth_min(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                         real_t smth_coef);
inline Real2D smooth_min(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                         real_t smth_coef);
inline Real2D smooth_min(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                         real_t smth_coef);
KOKKOS_INLINE_FUNCTION real_t smooth_min(real_t input_var1, real_t input_var2, real_t smth_coef);
inline Real2D smooth_max(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                         real_t smth_coef);
inline Real2D smooth_max(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                         real_t smth_coef);
inline Real2D smooth_max(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                         real_t smth_coef);
inline Real1D smooth_max(int ngrdcol, Real1D input_var1, real_t input_var2, real_t smth_coef);
KOKKOS_INLINE_FUNCTION real_t smooth_max(real_t input_var1, real_t input_var2, real_t smth_coef);

// =============================================================================
inline void calc_ddzt_umvm_sqd(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D um,
                               Real2D vm,           // In
                               Real2D ddzt_umvm_sqd // Out
) {
  // Description:
  // Computes the squared vertical norm of the derivative of the mean
  // horizontal wind.
  // References:
  // None
  // ---------------------------------------------------------------------------

  // --------------------------- Input Variables ---------------------------
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // um: Mean u wind on thermodynamic levels [m/s]
  // vm: Mean v wind on thermodynamic levels [m/s]
  // --------------------------- Output Variables --------------------------
  // ddzt_umvm_sqd: Squared vertical shear of horizontal wind [s^-2]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D ddzt_um("calc_ddzt_umvm_sqd::ddzt_um", ngrdcol,
                 nzm); // Vertical derivative of mean u wind [s^-1]
  Real2D ddzt_vm("calc_ddzt_umvm_sqd::ddzt_vm", ngrdcol,
                 nzm); // Vertical derivative of mean v wind [s^-1]

  Kokkos::deep_copy(ddzt_um, ddzt(nzm, nzt, ngrdcol, gr, um));
  Kokkos::deep_copy(ddzt_vm, ddzt(nzm, nzt, ngrdcol, gr, vm));
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_ddzt_umvm_sqd::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          ddzt_umvm_sqd(i - 1, k - 1) = ((ddzt_um(i - 1, k - 1) * ddzt_um(i - 1, k - 1)) +
                                         (ddzt_vm(i - 1, k - 1) * ddzt_vm(i - 1, k - 1)));
        });
  }
  return;
}

// =============================================================================
inline void calc_wp3_on_wp2(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D wp2,
                            Real2D wp3,                             // In
                            Real2D wp3_on_wp2, Real2D wp3_on_wp2_zt // Out
) {
  // Description:
  // Computes a smoothed ratio of w'^3 to w'^2 on thermodynamic and momentum
  // levels.
  // References:
  // None
  // ---------------------------------------------------------------------------

  // --------------------------- Input Variables ---------------------------
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // wp2: Variance of vertical velocity on momentum levels [m^2/s^2]
  // wp3: Third moment of vertical velocity on thermodynamic levels [m^3/s^3]
  // --------------------------- Output Variables --------------------------
  // wp3_on_wp2: Ratio of wp3 to wp2 on momentum levels [m/s]
  // wp3_on_wp2_zt: Ratio of wp3 to wp2 on thermodynamic levels [m/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D wp2_zt("calc_wp3_on_wp2::wp2_zt", ngrdcol,
                nzt); // Variance of vertical velocity on thermodynamic levels [m^2/s^2]

  Kokkos::deep_copy(wp2_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd));
  // Compute wp3 / wp2 on zt levels. Always use the interpolated value in the
  // denominator since it is less likely to create spikes.
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "calc_wp3_on_wp2::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          wp3_on_wp2_zt(i - 1, k - 1) =
              (wp3(i - 1, k - 1) /
               (wp2_zt(i - 1, k - 1) > w_tol_sqd ? wp2_zt(i - 1, k - 1) : w_tol_sqd));
          if ((wp3_on_wp2_zt(i - 1, k - 1) < real_t(0.))) {
            wp3_on_wp2_zt(i - 1, k - 1) =
                (-(real_t(1000.)) > wp3_on_wp2_zt(i - 1, k - 1) ? -(real_t(1000.))
                                                                : wp3_on_wp2_zt(i - 1, k - 1));
          } else {
            wp3_on_wp2_zt(i - 1, k - 1) =
                (real_t(1000.) < wp3_on_wp2_zt(i - 1, k - 1) ? real_t(1000.)
                                                             : wp3_on_wp2_zt(i - 1, k - 1));
          }
        });
  }
  Kokkos::deep_copy(wp3_on_wp2, zt2zm_api(nzm, nzt, ngrdcol, gr, wp3_on_wp2_zt));
  Kokkos::deep_copy(wp3_on_wp2_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, wp3_on_wp2));
  return;
}

// =============================================================================
inline void set_boundary_conditions_lhs(int diag_index, int low_bound, int high_bound, // In
                                        Real2D lhs,                                    // InOut
                                        std::optional<int> diag_index2,
                                        std::optional<int> low_bound2, // Optional In
                                        std::optional<int> high_bound2 // Optional In
) {
  // Set lower/upper boundary columns in a LAPACK band matrix. diag_index
  // selects the diagonal row. All input indices retain Native one-based values. Optional
  // second-variable boundaries require diag_index2. References: none.
  if ((low_bound2.has_value() || high_bound2.has_value()) && !diag_index2.has_value()) {
    throw std::runtime_error("Boundary index provided without diag_index.");
  }
  const bool has_low2 = low_bound2.has_value();
  const bool has_high2 = high_bound2.has_value();
  const int diag2 = diag_index2.value_or(1);
  const int low2 = low_bound2.value_or(1);
  const int high2 = high_bound2.value_or(1);
  // Each matrix row owns its boundary writes. Keeping all four stages in each
  // row preserves Native assignment order when boundary indices overlap.
  Kokkos::parallel_for(
      "set_boundary_conditions_lhs", lhs.extent_int(0), KOKKOS_LAMBDA(const int row) {
        // Set the lower boundaries for the first variable.
        lhs(row, low_bound - 1) = zero;
        if (row == diag_index - 1)
          lhs(row, low_bound - 1) = one;
        // Set the upper boundaries for the first variable.
        lhs(row, high_bound - 1) = zero;
        if (row == diag_index - 1)
          lhs(row, high_bound - 1) = one;
        // Set the lower boundaries for the second variable, if it is provided.
        if (has_low2) {
          lhs(row, low2 - 1) = zero;
          if (row == diag2 - 1)
            lhs(row, low2 - 1) = one;
        }
        // Set the upper boundaries for the second variable, if it is provided.
        if (has_high2) {
          lhs(row, high2 - 1) = zero;
          if (row == diag2 - 1)
            lhs(row, high2 - 1) = one;
        }
      });
}

// =============================================================================
inline void set_boundary_conditions_rhs(real_t low_value, int low_bound, real_t high_value,
                                        int high_bound, // In
                                        Real1D rhs,     // InOut
                                        std::optional<real_t> low_value2,
                                        std::optional<int> low_bound2, // Optional In
                                        std::optional<real_t> high_value2,
                                        std::optional<int> high_bound2 // Optional In
) {
  // Set the lower/upper LAPACK right-hand side values. Boundary indices retain
  // Native one-based values. Each optional bound requires its matching value.
  // References: none.
  // Stop execution if a boundary was provided without a value.
  if ((low_bound2.has_value() && !low_value2.has_value()) ||
      (high_bound2.has_value() && !high_value2.has_value())) {
    throw std::runtime_error("Boundary condition provided without value.");
  }
  const bool has_low2 = low_bound2.has_value();
  const bool has_high2 = high_bound2.has_value();
  const int low2 = low_bound2.value_or(1);
  const int high2 = high_bound2.value_or(1);
  const real_t value_low2 = low_value2.value_or(zero);
  const real_t value_high2 = high_value2.value_or(zero);
  // One work item retains lower1/upper1/lower2/upper2 overwrite order even
  // when the caller supplies overlapping boundary indices.
  Kokkos::parallel_for(
      "set_boundary_conditions_rhs", 1, KOKKOS_LAMBDA(const int) {
        // Set the lower and upper bounds for the first variable.
        rhs(low_bound - 1) = low_value;
        rhs(high_bound - 1) = high_value;
        // If a lower bound was given for the second variable, set it.
        if (has_low2)
          rhs(low2 - 1) = value_low2;
        // If an upper bound was given for the second variable, set it.
        if (has_high2)
          rhs(high2 - 1) = value_high2;
      });
}

// =============================================================================
inline void calc_stability_correction(int nzm, int ngrdcol, // In
                                      Real2D brunt_vaisala_freq_sqd, Real2D Lscale_zm,
                                      Real2D em,                     // In
                                      Real1D lambda0_stability_coef, // In
                                      Real2D stability_correction    // Out
) {
  // Description:
  // Compute one plus a bounded stable-stratification correction to the
  // turbulent length/time relation, using N^2, mixing length and TKE.
  // References:
  // --------------------------------------------------------------------

  // ---------------- Input Variables ----------------
  // brunt_vaisala_freq_sqd: Brunt-Vaisala frequency squared      [1/s^2]
  // Lscale_zm: Turbulent mixing length interp to m levs  [m]
  // em: Turbulent Kinetic Energy (TKE)            [m^2/s^2]
  // lambda0_stability_coef: CLUBB tunable parameter lambda0_stability_coef
  // ---------------- Output Variables ----------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.

  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_stability_correction::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          real_t lambda0_stability{};
          if ((brunt_vaisala_freq_sqd(i - 1, k - 1) > zero)) {
            lambda0_stability = lambda0_stability_coef(i - 1);
          } else {
            lambda0_stability = zero;
          }
          stability_correction(i - 1, k - 1) =
              (one + ((((lambda0_stability * brunt_vaisala_freq_sqd(i - 1, k - 1)) *
                        (Lscale_zm(i - 1, k - 1) * Lscale_zm(i - 1, k - 1))) /
                       em(i - 1, k - 1)) < three
                          ? (((lambda0_stability * brunt_vaisala_freq_sqd(i - 1, k - 1)) *
                              (Lscale_zm(i - 1, k - 1) * Lscale_zm(i - 1, k - 1))) /
                             em(i - 1, k - 1))
                          : three));
        });
  }
  return;
}

// =============================================================================
inline void calc_brunt_vaisala_freq_sqd(int nzm, int nzt, int ngrdcol, const grid& gr,
                                        Real2D thlm, // In
                                        Real2D exner, Real2D rtm, Real2D rcm, Real2D p_in_Pa,
                                        Real2D thvm,                          // In
                                        Real2D ice_supersat_frac,             // In
                                        int saturation_formula,               // In
                                        bool l_brunt_vaisala_freq_moist,      // In
                                        bool l_use_thvm_in_bv_freq,           // In
                                        bool l_modify_limiters_for_cnvg_test, // In
                                        Real1D bv_efold, real_t T0,           // In
                                        Real2D brunt_vaisala_freq_sqd,        // Out
                                        Real2D brunt_vaisala_freq_sqd_mixed,  // Out
                                        Real2D brunt_vaisala_freq_sqd_smth,   // Out
                                        f2k::stats_type* stats                // Optional InOut
) {
  // Description:
  // Calculate the Brunt-Vaisala frequency squared, N^2.
  // References:
  // ?
  // -----------------------------------------------------------------------

  // ---------------------------- Input Variables ----------------------------
  // thlm: th_l (thermo. levels)              [K]
  // exner: Exner function                     [-]
  // rtm: total water mixing ratio, r_t      [kg/kg]
  // rcm: cloud water mixing ratio, r_c      [kg/kg]
  // p_in_Pa: Air pressure                       [Pa]
  // thvm: Virtual potential temperature      [K]
  // saturation_formula: Integer that stores the saturation formula to be used
  // l_brunt_vaisala_freq_moist: Use a different formula for the Brunt-Vaisala frequency in
  // saturated atmospheres (from Durran and Klemp, 1982)
  // l_use_thvm_in_bv_freq: Use thvm in the calculation of Brunt-Vaisala frequency
  // Flag to activate modifications on limiters for convergence test
  // (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
  // (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
  // (reduce threshold on limiters for Ri_zm in mixing_length.F90)
  // bv_efold: Control parameter for inverse e-folding of
  // cloud fraction in the mixed Brunt Vaisala frequency
  // T0: Reference temperature (usually 300)  [K]
  // ---------------------------- Output Variables ----------------------------
  // brunt_vaisala_freq_sqd: Brunt-Vaisala frequency squared, N^2 [1/s^2]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D T_in_K("calc_brunt_vaisala_freq_sqd::T_in_K", ngrdcol, nzt);
  Real2D rsat("calc_brunt_vaisala_freq_sqd::rsat", ngrdcol, nzt);
  Real2D thm("calc_brunt_vaisala_freq_sqd::thm", ngrdcol, nzt);
  Real2D T_in_K_zm("calc_brunt_vaisala_freq_sqd::T_in_K_zm", ngrdcol, nzm);
  Real2D rsat_zm("calc_brunt_vaisala_freq_sqd::rsat_zm", ngrdcol, nzm);
  Real2D thm_zm("calc_brunt_vaisala_freq_sqd::thm_zm", ngrdcol, nzm);
  Real2D thvm_zm("calc_brunt_vaisala_freq_sqd::thvm_zm", ngrdcol, nzm);
  Real2D ice_supersat_frac_zm("calc_brunt_vaisala_freq_sqd::ice_supersat_frac_zm", ngrdcol, nzm);
  Real2D ddzt_thlm("calc_brunt_vaisala_freq_sqd::ddzt_thlm", ngrdcol, nzm);
  Real2D tmp_calc("calc_brunt_vaisala_freq_sqd::tmp_calc", ngrdcol, nzm);
  Real2D ddzt_thm("calc_brunt_vaisala_freq_sqd::ddzt_thm", ngrdcol, nzm);
  Real2D ddzt_rsat("calc_brunt_vaisala_freq_sqd::ddzt_rsat", ngrdcol, nzm);
  Real2D ddzt_rtm("calc_brunt_vaisala_freq_sqd::ddzt_rtm", ngrdcol, nzm);
  Real2D ddzt_thvm("calc_brunt_vaisala_freq_sqd::ddzt_thvm", ngrdcol, nzm);
  Real2D brunt_vaisala_freq_sqd_dry("calc_brunt_vaisala_freq_sqd::brunt_vaisala_freq_sqd_dry",
                                    ngrdcol, nzm);
  Real2D brunt_vaisala_freq_sqd_moist("calc_brunt_vaisala_freq_sqd::brunt_vaisala_freq_sqd_moist",
                                      ngrdcol, nzm);
  Real2D brunt_vaisala_freq_clipped("calc_brunt_vaisala_freq_sqd::brunt_vaisala_freq_clipped",
                                    ngrdcol, nzm);
  constexpr bool l_smooth_min_max = false;

  Kokkos::deep_copy(ddzt_thlm, ddzt(nzm, nzt, ngrdcol, gr, thlm));
  Kokkos::deep_copy(thvm_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, thvm, zero_threshold));
  Kokkos::deep_copy(ddzt_thvm, ddzt(nzm, nzt, ngrdcol, gr, thvm));
  // Dry Brunt-Vaisala frequency
  if (l_use_thvm_in_bv_freq) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_brunt_vaisala_freq_sqd::cells_0", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            brunt_vaisala_freq_sqd(i - 1, k - 1) =
                ((grav / thvm_zm(i - 1, k - 1)) * ddzt_thvm(i - 1, k - 1));
          });
    }
  } else {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_brunt_vaisala_freq_sqd::cells_1", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            brunt_vaisala_freq_sqd(i - 1, k - 1) = ((grav / T0) * ddzt_thlm(i - 1, k - 1));
          });
    }
  }
  Kokkos::deep_copy(T_in_K, thlm2T_in_K_api(nzt, ngrdcol, thlm, exner, rcm));
  Kokkos::deep_copy(T_in_K_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, T_in_K, zero_threshold));
  Kokkos::deep_copy(rsat,
                    sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula));
  Kokkos::deep_copy(rsat_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, rsat, zero_threshold));
  Kokkos::deep_copy(ddzt_rsat, ddzt(nzm, nzt, ngrdcol, gr, rsat));
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "calc_brunt_vaisala_freq_sqd::cells_2", ngrdcol * level_count,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          thm(i - 1, k - 1) =
              (thlm(i - 1, k - 1) + ((Lv / (Cp * exner(i - 1, k - 1))) * rcm(i - 1, k - 1)));
        });
  }
  Kokkos::deep_copy(thm_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, thm, zero_threshold));
  Kokkos::deep_copy(ddzt_thm, ddzt(nzm, nzt, ngrdcol, gr, thm));
  Kokkos::deep_copy(ddzt_rtm, ddzt(nzm, nzt, ngrdcol, gr, rtm));
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_brunt_vaisala_freq_sqd::cells_3", ngrdcol * level_count,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          brunt_vaisala_freq_sqd_dry(i - 1, k - 1) =
              ((grav / thm_zm(i - 1, k - 1)) * ddzt_thm(i - 1, k - 1));
        });
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_brunt_vaisala_freq_sqd::cells_4", ngrdcol * level_count,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          // In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and
          // Klemp (1982)
          brunt_vaisala_freq_sqd_moist(i - 1, k - 1) =
              (grav *
               ((((one + ((Lv * rsat_zm(i - 1, k - 1)) / (Rd * T_in_K_zm(i - 1, k - 1)))) /
                  (one + (((ep * (Lv * Lv)) * rsat_zm(i - 1, k - 1)) /
                          ((Cp * Rd) * (T_in_K_zm(i - 1, k - 1) * T_in_K_zm(i - 1, k - 1)))))) *
                 (((one / thm_zm(i - 1, k - 1)) * ddzt_thm(i - 1, k - 1)) +
                  ((Lv / (Cp * T_in_K_zm(i - 1, k - 1))) * ddzt_rsat(i - 1, k - 1)))) -
                ddzt_rtm(i - 1, k - 1)));
        });
  }
  // k=1, nzm
  Kokkos::deep_copy(ice_supersat_frac_zm,
                    zt2zm_api(nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold));
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_brunt_vaisala_freq_sqd::cells_5", ngrdcol * level_count,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          brunt_vaisala_freq_sqd_mixed(i - 1, k - 1) =
              (brunt_vaisala_freq_sqd_moist(i - 1, k - 1) +
               (Kokkos::exp((-(bv_efold(i - 1)) * ice_supersat_frac_zm(i - 1, k - 1))) *
                (brunt_vaisala_freq_sqd_dry(i - 1, k - 1) -
                 brunt_vaisala_freq_sqd_moist(i - 1, k - 1))));
        });
  }
  // The min function below smooths the slope discontinuity in brunt freq
  // and thereby allows tau to remain large in Sc layers in which thlm may
  // be slightly stably stratified.
  if (l_modify_limiters_for_cnvg_test) {
    // Remove the limiters to improve the solution convergence
    Kokkos::deep_copy(brunt_vaisala_freq_sqd_smth,
                      zm2zt2zm(nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_sqd_mixed));
  } else {
    // default method
    if (l_smooth_min_max) {
      {
        const int level_count = nzm;
        Kokkos::parallel_for(
            "calc_brunt_vaisala_freq_sqd::cells_6", ngrdcol * level_count,
            KOKKOS_LAMBDA(const int cell) {
              const int i = cell % ngrdcol + 1;
              const int k = cell / ngrdcol + 1;
              tmp_calc(i - 1, k - 1) =
                  (real_t(1.e8) * ((Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1)) *
                                    Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1))) *
                                   Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1))));
            });
      }
      Kokkos::deep_copy(brunt_vaisala_freq_clipped,
                        smooth_min(nzm, ngrdcol, brunt_vaisala_freq_sqd_mixed, tmp_calc,
                                   (real_t(1.0e-4) * min_max_smth_mag)));
      Kokkos::deep_copy(brunt_vaisala_freq_sqd_smth,
                        zm2zt2zm(nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_clipped));
    } else {
      // l_smooth_min_max
      {
        const int level_count = nzm;
        Kokkos::parallel_for(
            "calc_brunt_vaisala_freq_sqd::cells_7", ngrdcol * level_count,
            KOKKOS_LAMBDA(const int cell) {
              const int i = cell % ngrdcol + 1;
              const int k = cell / ngrdcol + 1;
              brunt_vaisala_freq_clipped(i - 1, k - 1) =
                  (brunt_vaisala_freq_sqd_mixed(i - 1, k - 1) <
                           (real_t(1.e8) *
                            ((Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1)) *
                              Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1))) *
                             Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1))))
                       ? brunt_vaisala_freq_sqd_mixed(i - 1, k - 1)
                       : (real_t(1.e8) *
                          ((Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1)) *
                            Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1))) *
                           Kokkos::abs(brunt_vaisala_freq_sqd_mixed(i - 1, k - 1)))));
            });
      }
      Kokkos::deep_copy(brunt_vaisala_freq_sqd_smth,
                        zm2zt2zm(nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_clipped));
    }
  }
  if (l_brunt_vaisala_freq_moist) {
    Kokkos::deep_copy(brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_moist);
  }
  if (stats != nullptr) {
    if (stats->flags().l_sample) {
      f2k::stats_update("bv_freq_sqd", brunt_vaisala_freq_sqd, // In
                        *stats                                 // InOut
      );
      f2k::stats_update("bv_freq_sqd_mixed", brunt_vaisala_freq_sqd_mixed, // In
                        *stats                                             // InOut
      );
      f2k::stats_update("bv_freq_sqd_moist", brunt_vaisala_freq_sqd_moist, // In
                        *stats                                             // InOut
      );
      f2k::stats_update("bv_freq_sqd_dry", brunt_vaisala_freq_sqd_dry, // In
                        *stats                                         // InOut
      );
      f2k::stats_update("bv_freq_sqd_smth", brunt_vaisala_freq_sqd_smth, // In
                        *stats                                           // InOut
      );
    }
  }
  return;
}

// =============================================================================
inline void compute_Cx_fnc_Richardson(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                                      Real2D Lscale_zm, Real2D ddzt_umvm_sqd,
                                      Real2D rho_ds_zm,                     // In
                                      Real2D brunt_vaisala_freq_sqd,        // In
                                      Real2D brunt_vaisala_freq_sqd_mixed,  // In
                                      Real2D clubb_params,                  // In
                                      bool l_use_shear_Richardson,          // In
                                      bool l_modify_limiters_for_cnvg_test, // In
                                      Real2D Cx_fnc_Richardson              // Out
) {
  // Description:
  // Compute Cx as a function of the Richardson number
  // References:
  // cam:ticket:59
  // -----------------------------------------------------------------------

  // ------------------------------ Constant Parameters ------------------------------
  // ------------------------------ Input Variables ------------------------------
  // brunt_vaisala_freq_sqd: Brunt-Vaisala frequency squared, N^2  [1/s^2]
  // brunt_vaisala_freq_sqd_mixed: A mixture of dry and moist N^2        [s^-2]
  // ddzt_umvm_sqd: Squared vertical norm of derivative of
  // mean horizontal wind speed            [s^-2]
  // Lscale_zm: Length scale on momentum levels       [m]
  // rho_ds_zm: Dry static density on momentum levels [kg/m^3]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // l_use_shear_Richardson: Use shear in the calculation of Richardson number
  // Flag to activate modifications on limiters for convergence test
  // (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
  // (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
  // (reduce threshold on limiters for sqrt_Ri_zm in mixing_length.F90)
  // ------------------------------ Output Variable ------------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  constexpr real_t Richardson_num_divisor_threshold = real_t(1.0e-6);
  constexpr real_t Cx_fnc_Richardson_below_ground_value = one;
  constexpr bool l_Cx_fnc_Richardson_vert_avg =
      false; // Vertically average Cx_fnc_Richardson over a distance of Lscale.
  Real2D Ri_zm_Cx("compute_Cx_fnc_Richardson::Ri_zm_Cx", ngrdcol, nzm);
  Real2D fnc_Richardson("compute_Cx_fnc_Richardson::fnc_Richardson", ngrdcol, nzm);
  Real2D fnc_Richardson_clipped("compute_Cx_fnc_Richardson::fnc_Richardson_clipped", ngrdcol, nzm);
  Real2D fnc_Richardson_smooth("compute_Cx_fnc_Richardson::fnc_Richardson_smooth", ngrdcol, nzm);
  Real2D Cx_fnc_interp("compute_Cx_fnc_Richardson::Cx_fnc_interp", ngrdcol, nzm);
  Real2D Cx_fnc_Richardson_avg("compute_Cx_fnc_Richardson::Cx_fnc_Richardson_avg", ngrdcol, nzm);
  real_t invrs_num_div_thresh{};
  int smth_type = 1;

  invrs_num_div_thresh = (one / Richardson_num_divisor_threshold);
  if (l_use_shear_Richardson) {
    calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_sqd_mixed, ddzt_umvm_sqd, // In
               real_t(1.0e-7), real_t(1.0e-7),                            // In
               Ri_zm_Cx                                                   // Out
    );
  } else {
    // Note1: We kind of want this calculation to be done in calc_Ri_zm, as well.
    // But multiplying by the inverse and dividing are not numerically equal.
    // So to preserve BFBness, we decided to keep this as is, for now.
    // Note2: Keep in mind, that this Brunt Vaisala variable can contain negative values.
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "compute_Cx_fnc_Richardson::cells_0", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            Ri_zm_Cx(i - 1, k - 1) = (brunt_vaisala_freq_sqd(i - 1, k - 1) * invrs_num_div_thresh);
          });
    }
  }
  // Cx_fnc_Richardson is interpolated based on the value of Richardson_num
  // The min function ensures that Cx does not exceed Cx_max, regardless of the
  // value of Richardson_num_max.
  if (l_modify_limiters_for_cnvg_test) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "compute_Cx_fnc_Richardson::cells_1", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            real_t invrs_min_max_diff{};
            real_t Richardson_num_max{};
            real_t Richardson_num_min{};
            Richardson_num_max = clubb_params(i - 1, iRichardson_num_max);
            Richardson_num_min = clubb_params(i - 1, iRichardson_num_min);
            invrs_min_max_diff = (one / (Richardson_num_max - Richardson_num_min));
            fnc_Richardson(i - 1, k - 1) =
                ((Ri_zm_Cx(i - 1, k - 1) - clubb_params(i - 1, iRichardson_num_min)) *
                 invrs_min_max_diff);
          });
    }
    Kokkos::deep_copy(fnc_Richardson_clipped,
                      smooth_min(nzm, ngrdcol, one, fnc_Richardson, min_max_smth_mag));
    Kokkos::deep_copy(fnc_Richardson_smooth,
                      smooth_max(nzm, ngrdcol, zero, fnc_Richardson_clipped, min_max_smth_mag));
    // use smoothed max amd min to achive smoothed profile and avoid discontinuities
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "compute_Cx_fnc_Richardson::cells_2", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            real_t Cx_max{};
            real_t Cx_min{};
            Cx_max = clubb_params(i - 1, iCx_max);
            Cx_min = clubb_params(i - 1, iCx_min);
            Cx_fnc_interp(i - 1, k - 1) =
                ((fnc_Richardson_smooth(i - 1, k - 1) * (Cx_max - Cx_min)) + Cx_min);
          });
    }
    Kokkos::deep_copy(Cx_fnc_Richardson, zm2zt2zm(nzm, nzt, ngrdcol, gr, Cx_fnc_interp));
  } else {
    // default method
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "compute_Cx_fnc_Richardson::cells_3", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            real_t invrs_min_max_diff{};
            real_t Richardson_num_max{};
            real_t Richardson_num_min{};
            real_t Cx_max{};
            real_t Cx_min{};
            Richardson_num_max = clubb_params(i - 1, iRichardson_num_max);
            Richardson_num_min = clubb_params(i - 1, iRichardson_num_min);
            Cx_max = clubb_params(i - 1, iCx_max);
            Cx_min = clubb_params(i - 1, iCx_min);
            invrs_min_max_diff = (one / (Richardson_num_max - Richardson_num_min));
            Cx_fnc_Richardson(i - 1, k - 1) =
                ((((((Richardson_num_max < Ri_zm_Cx(i - 1, k - 1)
                          ? Richardson_num_max
                          : Ri_zm_Cx(i - 1, k - 1)) > Richardson_num_min
                         ? (Richardson_num_max < Ri_zm_Cx(i - 1, k - 1) ? Richardson_num_max
                                                                        : Ri_zm_Cx(i - 1, k - 1))
                         : Richardson_num_min) -
                    Richardson_num_min) *
                   invrs_min_max_diff) *
                  (Cx_max - Cx_min)) +
                 Cx_min);
          });
    }
  }
  // Native keeps this parameter false. Its dormant block subsequently copies
  // uninitialized Cx_fnc_Richardson_avg over the computed result; retained here
  // for source fidelity rather than silently changing that Native branch.
  if (l_Cx_fnc_Richardson_vert_avg) {
    Kokkos::deep_copy(Cx_fnc_Richardson,
                      Lscale_width_vert_avg(nzm, ngrdcol, gr, smth_type, Cx_fnc_Richardson,
                                            Lscale_zm, rho_ds_zm,
                                            Cx_fnc_Richardson_below_ground_value));
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "compute_Cx_fnc_Richardson::cells_4", ngrdcol * level_count,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            Cx_fnc_Richardson(i - 1, k - 1) = Cx_fnc_Richardson_avg(i - 1, k - 1);
          });
    }
  }
  // On some compilers, roundoff error can result in Cx_fnc_Richardson being
  // slightly outside the range [0,1]. Thus, it is clipped here.
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "compute_Cx_fnc_Richardson::cells_5", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          Cx_fnc_Richardson(i - 1, k - 1) =
              (zero > (one < Cx_fnc_Richardson(i - 1, k - 1) ? one
                                                             : Cx_fnc_Richardson(i - 1, k - 1))
                   ? zero
                   : (one < Cx_fnc_Richardson(i - 1, k - 1) ? one
                                                            : Cx_fnc_Richardson(i - 1, k - 1)));
        });
  }
  return;
}

// =============================================================================
inline void calc_Ri_zm(int nzm, int ngrdcol, Real2D bv_freq_sqd, Real2D shear, real_t lim_bv, // In
                       real_t lim_shear,                                                      // In
                       Real2D Ri_zm                                                           // Out
) {
  // Description:
  // Calculate the Richardson number as a quotient of a clipped Brunt-Vaisala frequency
  // and a clipped shear
  // References:
  // ?

  // -------------------------- Input Variables --------------------------
  // bv_freq_sqd: Square of Brunt Vaisala frequency
  // shear: Shear variable, usually norm squared of horizontal wind speeds
  // lim_bv: Clipping value in numerator
  // lim_shear: Clipping value in denominator
  // -------------------------- Result Variable --------------------------
  // Ri_zm: Richardson number

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.

  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_Ri_zm::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          Ri_zm(i - 1, k - 1) =
              ((bv_freq_sqd(i - 1, k - 1) > lim_bv ? bv_freq_sqd(i - 1, k - 1) : lim_bv) /
               (shear(i - 1, k - 1) > lim_shear ? shear(i - 1, k - 1) : lim_shear));
        });
  }
}

// =============================================================================
inline Real2D Lscale_width_vert_avg(int nzm, int ngrdcol, const grid& gr, int smth_type,    // In
                                    Real2D var_profile, Real2D Lscale_zm, Real2D rho_ds_zm, // In
                                    real_t var_below_ground_value                           // In
) {
  // Description:
  // Averages a profile with a running mean of width Lscale_zm
  // References:
  // cam:ticket:59

  // -------------------------- Input Variables --------------------------
  // var_profile: Profile on momentum levels
  // Lscale_zm: Lscale on momentum levels
  // rho_ds_zm: Dry static density on momentum levels [kg/m^3]; the Native
  // declaration labels it energy, but the weighted integral uses density.
  // var_below_ground_value: Value to use below ground
  // -------------------------- Result Variable --------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D Lscale_width_vert_avg_output("Lscale_width_vert_avg::Lscale_width_vert_avg_output",
                                      ngrdcol,
                                      nzm); // Vertically averaged profile (on momentum levels)
  Real2D one_half_avg_width("Lscale_width_vert_avg::one_half_avg_width", ngrdcol, nzm);
  Real2D numer_terms("Lscale_width_vert_avg::numer_terms", ngrdcol, nzm);
  Real2D denom_terms("Lscale_width_vert_avg::denom_terms", ngrdcol, nzm);
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_grid_dir = gr.grid_dir;
  const auto gr_dzm = gr.dzm;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_zm = gr.zm;

  if ((smth_type == 1)) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "Lscale_width_vert_avg::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            one_half_avg_width(i - 1, k - 1) =
                (Lscale_zm(i - 1, k - 1) > real_t(500.0) ? Lscale_zm(i - 1, k - 1) : real_t(500.0));
          });
    }
  } else if ((smth_type == 2)) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "Lscale_width_vert_avg::cells_1", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            one_half_avg_width(i - 1, k - 1) = real_t(60.0);
          });
    }
  }
  // Pre calculate numerator and denominator terms
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "Lscale_width_vert_avg::cells_2", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          numer_terms(i - 1, k - 1) =
              (((rho_ds_zm(i - 1, k - 1) * gr_grid_dir) * gr_dzm(i - 1, k - 1)) *
               var_profile(i - 1, k - 1));
          denom_terms(i - 1, k - 1) =
              ((rho_ds_zm(i - 1, k - 1) * gr_grid_dir) * gr_dzm(i - 1, k - 1));
        });
  }
  // For every grid level
  {
    const int level_count =
        ((gr_grid_dir_indx) > 0 ? (gr_k_lb_zm) <= (gr_k_ub_zm) : (gr_k_lb_zm) >= (gr_k_ub_zm))
            ? ((gr_k_ub_zm) - (gr_k_lb_zm)) / (gr_grid_dir_indx) + 1
            : 0;
    Kokkos::parallel_for(
        "Lscale_width_vert_avg::cells_3", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = (gr_k_lb_zm) + (cell / ngrdcol) * (gr_grid_dir_indx);
          int k_avg_lower{};
          int k_avg_upper{};
          int k_avg{};
          int n_below_ground_levels{};
          real_t numer_integral{};
          real_t denom_integral{};
          // -----------------------------------------------------------------------
          // Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
          // Note: Outdated explanation of version that improves CPU performance
          // below. Reworked due to it requiring a k dependency. Now we
          // begin looking for k_avg_upper and k_avg_lower starting at
          // the kth level.
          // Outdated but potentially useful note:
          // k_avg_upper and k_avg_lower can be saved each loop iteration, this
          // reduces iterations beacuse the kth values are likely to be within
          // one or two grid levels of the k-1th values. Less searching is required
          // by starting the search at the previous values and incrementing or
          // decrement as needed.
          // -----------------------------------------------------------------------
          // Determine if k_avg_upper needs to increment
          for (k_avg_upper = k;
               ((gr_grid_dir_indx) > 0 ? k_avg_upper <= ((gr_k_ub_zm - gr_grid_dir_indx))
                                       : k_avg_upper >= ((gr_k_ub_zm - gr_grid_dir_indx)));
               k_avg_upper += gr_grid_dir_indx) {
            if (((gr_zm(i - 1, (k_avg_upper + gr_grid_dir_indx) - 1) - gr_zm(i - 1, k - 1)) >
                 one_half_avg_width(i - 1, k - 1))) {
              break;
            }
          }
          // Determine if k_avg_lower needs to decrement
          for (k_avg_lower = k;
               ((-(gr_grid_dir_indx)) > 0 ? k_avg_lower <= ((gr_k_lb_zm + gr_grid_dir_indx))
                                          : k_avg_lower >= ((gr_k_lb_zm + gr_grid_dir_indx)));
               k_avg_lower += -(gr_grid_dir_indx)) {
            if (((gr_zm(i - 1, k - 1) - gr_zm(i - 1, (k_avg_lower - gr_grid_dir_indx) - 1)) >
                 one_half_avg_width(i - 1, k - 1))) {
              break;
            }
          }
          // Compute the number of levels below ground to include.
          if (((gr_grid_dir_indx * k_avg_lower) > (gr_grid_dir_indx * gr_k_lb_zm))) {
            // k=gr%k_lb_zm, the lowest level, is not included in the average, so no
            // below-ground levels should be included.
            n_below_ground_levels = 0;
            numer_integral = zero;
            denom_integral = zero;
          } else {
            // The number of below-ground levels included is equal to the distance
            // below the lowest level spanned by one_half_avg_width(k)
            // divided by the distance between vertical levels below ground; the
            // latter is assumed to be the same as the distance between the first and
            // second vertical levels.
            n_below_ground_levels =
                static_cast<int>(((one_half_avg_width(i - 1, k - 1) -
                                   (gr_zm(i - 1, k - 1) - gr_zm(i - 1, gr_k_lb_zm - 1))) /
                                  (gr_zm(i - 1, (gr_k_lb_zm + gr_grid_dir_indx) - 1) -
                                   gr_zm(i - 1, gr_k_lb_zm - 1))));
            numer_integral = ((n_below_ground_levels * denom_terms(i - 1, gr_k_lb_zm - 1)) *
                              var_below_ground_value);
            denom_integral = (n_below_ground_levels * denom_terms(i - 1, gr_k_lb_zm - 1));
          }
          // Add numerator and denominator terms for all above-ground levels
          for (k_avg = k_avg_lower;
               ((gr_grid_dir_indx) > 0 ? k_avg <= (k_avg_upper) : k_avg >= (k_avg_upper));
               k_avg += gr_grid_dir_indx) {
            numer_integral = (numer_integral + numer_terms(i - 1, k_avg - 1));
            denom_integral = (denom_integral + denom_terms(i - 1, k_avg - 1));
          }
          Lscale_width_vert_avg_output(i - 1, k - 1) = (numer_integral / denom_integral);
        });
  }
  return Lscale_width_vert_avg_output;
}

// =============================================================================
inline void wp23_term_splat_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,
                                Real1D C_wp2_splat, // In
                                Real2D brunt_vaisala_freq_sqd_mixed, Real2D Lscale_zm,
                                Real2D rho_ds_zm,                          // In
                                f2k::stats_type& stats,                    // InOut
                                Real2D lhs_splat_wp2, Real2D lhs_splat_wp3 // Out
) {
  // Description:
  // Calculate the LHS coefficients for the wp2 and wp3 splatting terms.
  // References:
  // -----------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // brunt_vaisala_freq_sqd_mixed: Mixture of dry and moist Brunt-Vaisala frequency squared [s^-2]
  // Lscale_zm: Turbulent mixing length on momentum levels [m]
  // rho_ds_zm: Dry, static density on momentum levels [kg/m^3]
  // C_wp2_splat: Model parameter C_wp2_splat             [ -]
  // --------------------- Input/Output Variables ---------------------
  // stats: Statistics accumulator
  // --------------------- Output Variable ---------------------
  // lhs_splat_wp2: LHS coefficient of wp2 splatting term  [1/s]
  // lhs_splat_wp3: LHS coefficient of wp3 splatting term [1/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D brunt_vaisala_freq_sqd_splat(
      "wp23_term_splat_lhs::brunt_vaisala_freq_sqd_splat", ngrdcol,
      nzm); // Vertically averaged Brunt-Vaisala frequency squared [s^-2]
  Real2D brunt_vaisala_freq_splat_clipped("wp23_term_splat_lhs::brunt_vaisala_freq_splat_clipped",
                                          ngrdcol,
                                          nzm); // Nonnegative Brunt-Vaisala frequency [1/s]
  Real2D brunt_vaisala_freq_splat_smooth("wp23_term_splat_lhs::brunt_vaisala_freq_splat_smooth",
                                         ngrdcol,
                                         nzm); // Smoothed nonnegative Brunt-Vaisala frequency [1/s]
  Real2D brunt_vaisala_freq_splat_clipped_zt(
      "wp23_term_splat_lhs::brunt_vaisala_freq_splat_clipped_zt", ngrdcol,
      nzt); // Nonnegative Brunt-Vaisala frequency on thermo. levels [1/s]
  constexpr real_t below_grnd_val =
      real_t(0.01);            // Below-ground value for vertical averaging [s^-2]
  constexpr int smth_type = 2; // Lscale_width_vert_avg smoothing type

  // Calculate Brunt-Vaisala frequency used for splatting.
  Kokkos::deep_copy(brunt_vaisala_freq_sqd_splat,
                    Lscale_width_vert_avg(nzm, ngrdcol, gr, smth_type, brunt_vaisala_freq_sqd_mixed,
                                          Lscale_zm, rho_ds_zm, below_grnd_val));
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "wp23_term_splat_lhs::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          brunt_vaisala_freq_splat_clipped(i - 1, k - 1) =
              Kokkos::sqrt((zero > brunt_vaisala_freq_sqd_splat(i - 1, k - 1)
                                ? zero
                                : brunt_vaisala_freq_sqd_splat(i - 1, k - 1)));
        });
  }
  Kokkos::deep_copy(brunt_vaisala_freq_splat_smooth,
                    zm2zt2zm(nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_splat_clipped));
  Kokkos::deep_copy(brunt_vaisala_freq_splat_clipped_zt,
                    zm2zt_api(nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_splat_clipped));
  // Vertical compression of eddies causes gustiness (increase in up2 and vp2).
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "wp23_term_splat_lhs::cells_1", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          lhs_splat_wp2(i - 1, k - 1) =
              (+(C_wp2_splat(i - 1)) * brunt_vaisala_freq_splat_smooth(i - 1, k - 1));
        });
  }
  // Vertical compression of eddies also diminishes w'3.
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "wp23_term_splat_lhs::cells_2", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          lhs_splat_wp3(i - 1, k - 1) = (((+(one_half)*three) * C_wp2_splat(i - 1)) *
                                         brunt_vaisala_freq_splat_clipped_zt(i - 1, k - 1));
        });
  }
  if (stats.flags().l_sample) {
    f2k::stats_update("bv_freq_sqd_splat", brunt_vaisala_freq_sqd_splat, // In
                      stats                                              // InOut
    );
  }
  return;
}

// =============================================================================
inline Real2D smooth_min_sclr_idx(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                                  real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the min function, using one scalar and
  // one column/level array as inputs. For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // input_var2: Units vary
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_min_sclr_idx::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_min_sclr_idx::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1 + input_var2(i - 1, k - 1)) -
                           Kokkos::sqrt((((input_var1 - input_var2(i - 1, k - 1)) *
                                          (input_var1 - input_var2(i - 1, k - 1))) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
inline Real2D smooth_min_array_scalar(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                                      real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the min function, using one scalar and
  // one column/level array as inputs. For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_min_array_scalar::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_min_array_scalar::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1(i - 1, k - 1) + input_var2) -
                           Kokkos::sqrt((((input_var1(i - 1, k - 1) - input_var2) *
                                          (input_var1(i - 1, k - 1) - input_var2)) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
inline Real2D smooth_min_arrays(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                                real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed minimum of two column/level arrays.
  // For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables-----------------------------
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_min_arrays::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_min_arrays::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1(i - 1, k - 1) + input_var2(i - 1, k - 1)) -
                           Kokkos::sqrt((((input_var1(i - 1, k - 1) - input_var2(i - 1, k - 1)) *
                                          (input_var1(i - 1, k - 1) - input_var2(i - 1, k - 1))) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t smooth_min_scalars(real_t input_var1, real_t input_var2,
                                                 real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the min function, using two scalars as inputs.
  // For more details, see the interface in this file.
  // References:
  // See clubb:ticket: 965
  // ----------------------------------------------------------------------

  // Input Variables
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // Output Variables

  real_t output_var{}; // Same unit as input_var1 and input_var2

  output_var = (one_half * ((input_var1 + input_var2) -
                            Kokkos::sqrt((((input_var1 - input_var2) * (input_var1 - input_var2)) +
                                          (smth_coef * smth_coef)))));
  return output_var;
}

// =============================================================================
inline Real2D smooth_max_sclr_idx(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                                  real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the max function, using one scalar and
  // one column/level array as inputs. For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // input_var2: Units vary
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_max_sclr_idx::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_max_sclr_idx::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1 + input_var2(i - 1, k - 1)) +
                           Kokkos::sqrt((((input_var1 - input_var2(i - 1, k - 1)) *
                                          (input_var1 - input_var2(i - 1, k - 1))) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
inline Real2D smooth_max_array_scalar(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                                      real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the max function, using one scalar and
  // one column/level array as inputs. For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_max_array_scalar::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_max_array_scalar::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1(i - 1, k - 1) + input_var2) +
                           Kokkos::sqrt((((input_var1(i - 1, k - 1) - input_var2) *
                                          (input_var1(i - 1, k - 1) - input_var2)) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
inline Real1D smooth_max_array_1D_scalar(int ngrdcol, Real1D input_var1, real_t input_var2,
                                         real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed version of the max function, using one scalar and
  // one column array as inputs. For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real1D output_var("smooth_max_array_1D_scalar::output_var",
                    ngrdcol); // Same unit as input_var1 and input_var2

  Kokkos::parallel_for(
      "smooth_max_array_1D_scalar::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column + 1;
        output_var(i - 1) =
            (one_half *
             ((input_var1(i - 1) + input_var2) +
              Kokkos::sqrt((((input_var1(i - 1) - input_var2) * (input_var1(i - 1) - input_var2)) +
                            (smth_coef * smth_coef)))));
      });
  return output_var;
}

// =============================================================================
inline Real2D smooth_max_arrays(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                                real_t smth_coef // In
) {
  // Description:
  // Computes a smoothed maximum of two column/level arrays.
  // For more details, see the interface in this file.
  // References:
  // See clubb:ticket:894, updated version: 965
  // ----------------------------------------------------------------------

  // ----------------------------- Input Variables -----------------------------
  // input_var1: Units vary
  // input_var2: Units vary
  // smth_coef: "intensity" of the smoothing. Should be of a similar magnitude to
  // that of the data structures input_var1 and input_var2
  // ----------------------------- Output Variables -----------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D output_var("smooth_max_arrays::output_var", ngrdcol,
                    nz); // Same unit as input_var1 and input_var2

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_max_arrays::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          output_var(i - 1, k - 1) =
              (one_half * ((input_var1(i - 1, k - 1) + input_var2(i - 1, k - 1)) +
                           Kokkos::sqrt((((input_var1(i - 1, k - 1) - input_var2(i - 1, k - 1)) *
                                          (input_var1(i - 1, k - 1) - input_var2(i - 1, k - 1))) +
                                         (smth_coef * smth_coef)))));
        });
  }
  return output_var;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t smooth_max_scalars(real_t input_var1, real_t input_var2,
                                                 real_t smth_coef // In
) {
  // Compute a smoothed maximum of two scalar values, returning the same units
  // as the inputs. smth_coef controls the smoothing and should have a similar
  // magnitude to input_var1 and input_var2. Reference: CLUBB ticket 965.
  real_t output_var; // Smoothed maximum [input units]

  //----------------------------- Begin Code -----------------------------
  output_var =
      one_half *
      ((input_var1 + input_var2) +
       Kokkos::sqrt((input_var1 - input_var2) * (input_var1 - input_var2) + smth_coef * smth_coef));
  return output_var;
}

// =============================================================================
inline Real2D smooth_heaviside_peskin(int nz, int ngrdcol, Real2D input, // In
                                      real_t smth_range                  // In
) {
  // Description:
  // Computes a smoothed heaviside function as in
  // [Lin, Lee et al., 2005, A level set characteristic Galerkin
  // finite element method for free surface flows], equation (2)
  // References:
  // See clubb:ticket:965
  // ----------------------------------------------------------------------

  // ------------------------- Input Variables -------------------------
  // input: Units vary
  // smth_range: Smooth Heaviside function on [-smth_range, smth_range]
  // ------------------------- Output Variables -------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D smth_output("smooth_heaviside_peskin::smth_output", ngrdcol, nz); // Same units as input

  {
    const int level_count = nz;
    Kokkos::parallel_for(
        "smooth_heaviside_peskin::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          real_t input_over_smth_range{};
          if ((input(i - 1, k - 1) < -(smth_range))) {
            smth_output(i - 1, k - 1) = zero;
          } else if ((input(i - 1, k - 1) > smth_range)) {
            smth_output(i - 1, k - 1) = one;
          } else {
            // Note that this case will only ever be reached if smth_range != 0,
            // so this division is fine and should not cause any issues
            input_over_smth_range = (input(i - 1, k - 1) / smth_range);
            smth_output(i - 1, k - 1) =
                (one_half * ((one + input_over_smth_range) +
                             (invrs_pi * Kokkos::sin((pi * input_over_smth_range)))));
          }
        });
  }
  return smth_output;
}

// =============================================================================
inline void calc_xpwp_1D(const grid& gr, Real1D Km_zm, Real1D xm, // In
                         Real1D xpwp                              // Out
) {
  // Description:
  // Compute x'w' from x<k>, x<k-1>, Km_zm and invrs_dzm
  // References:
  // None
  // -----------------------------------------------------------------------

  // ----------------------- Input variables -----------------------
  // Km_zm: Eddy diff. (k momentum level)                 [m^2/s]
  // xm: x (k thermo level)                            [units vary]
  // ----------------------- Output variable -----------------------
  // xpwp: x'w'   [(units vary)(m/s)]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_nzm = gr.nzm;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // Solve for x'w' at all intermediate model levels. Native DO has no
  // iterations when there are fewer than three momentum levels.
  const int interior_count = gr_nzm > 2 ? gr_nzm - 2 : 0;
  Kokkos::parallel_for(
      "calc_xpwp_1D::interior", interior_count, KOKKOS_LAMBDA(const int level) {
        const int k = level + 2;
        xpwp(k - 1) = ((Km_zm(k - 1) * gr_invrs_dzm(1 - 1, k - 1)) * (xm(k - 1) - xm((k - 1) - 1)));
      });
  return;
}

// =============================================================================
inline void calc_xpwp_2D(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                         Real2D Km_zm, Real2D xm,                       // In
                         Real2D xpwp                                    // Out
) {
  // Description:
  // Compute x'w' from x<k>, x<k-1>, Km_zm and invrs_dzm
  // References:
  // None
  // -----------------------------------------------------------------------

  // ----------------------- Input variables -----------------------
  // Km_zm: Eddy diff. (k momentum level)                 [m^2/s]
  // xm: x (k thermo level)                            [units vary]
  // ----------------------- Output variable -----------------------
  // xpwp: x'w'   [(units vary)(m/s)]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // Solve for x'w' at all intermediate model levels.
  {
    const int level_count =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "calc_xpwp_2D::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = (2) + (cell / ngrdcol) * (1);
          xpwp(i - 1, k - 1) = ((Km_zm(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1)) *
                                (xm(i - 1, k - 1) - xm(i - 1, (k - 1) - 1)));
        });
  }
  return;
}

// =============================================================================
template <class Rho, class Field, class Dz>
inline real_t vertical_avg(int total_idx, Rho rho_ds, Field field, Dz dz // In
) {
  // Description:
  // Computes the density-weighted vertical average of a field.
  //
  // The average value of a function, f, over a set domain, [a,b], is
  // calculated by the equation:
  //
  // f_avg = ( INT(a:b) f*g ) / ( INT(a:b) g );
  //
  // as long as f is continous and g is nonnegative and integrable.  Therefore,
  // the density-weighted (by dry, static, base-static density) vertical
  // average value of any model field, x, is calculated by the equation:
  //
  // x_avg|_z = ( INT(z_bot:z_top) x rho_ds dz )
  //            / ( INT(z_bot:z_top) rho_ds dz );
  //
  // where z_bot is the bottom of the vertical domain, and z_top is the top of
  // the vertical domain.
  //
  // This calculation is done slightly differently depending on whether x is a
  // thermodynamic-level or a momentum-level variable.
  //
  // Thermodynamic-level computation:
  //
  // For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
  // numerator integral, is calculated as:
  //
  // SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
  //
  // where k is the index of the given thermodynamic level, x and rho_ds are
  // both thermodynamic-level variables, and delta_z(k) = zm(k) - zm(k-1).  The
  // indices k_bot and k_top are the indices of the respective lower and upper
  // thermodynamic levels involved in the integration.
  //
  // Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
  // is calculated as:
  //
  // SUM(k_bot:k_top) rho_ds(k) delta_z(k).
  //
  // The first (k=1) thermodynamic level is below ground (or below the
  // official lower boundary at the first momentum level), so it should not
  // count in a vertical average, whether that vertical average is used for
  // the hole-filling scheme or for statistical purposes. Begin no lower
  // than level k=2, which is the first thermodynamic level above ground (or
  // above the model lower boundary).
  //
  // For cases where hole-filling over the entire (global) vertical domain
  // is desired, or where statistics over the entire (global) vertical
  // domain are desired, the lower (thermodynamic-level) index of k = 2 and
  // the upper (thermodynamic-level) index of k = gr%nz, means that the
  // overall vertical domain will be gr%zm(1,gr%nz) - gr%zm(1,1).
  //
  //
  // Momentum-level computation:
  //
  // For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
  // numerator integral, is calculated as:
  //
  // SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
  //
  // where k is the index of the given momentum level, x and rho_ds are both
  // momentum-level variables, and delta_z(k) = zt(k+1) - zt(k).  The indices
  // k_bot and k_top are the indices of the respective lower and upper momentum
  // levels involved in the integration.
  //
  // Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
  // is calculated as:
  //
  // SUM(k_bot:k_top) rho_ds(k) delta_z(k).
  //
  // The first (k=1) momentum level is right at ground level (or right at
  // the official lower boundary).  The momentum level variables that call
  // the hole-filling scheme have set values at the surface (or lower
  // boundary), and those set values should not be changed.  Therefore, the
  // vertical average (for purposes of hole-filling) should not include the
  // surface level (or lower boundary level).  For hole-filling purposes,
  // begin no lower than level k=2, which is the second momentum level above
  // ground (or above the model lower boundary).  Likewise, the value at the
  // model upper boundary (k=gr%nz) is also set for momentum level
  // variables.  That value should also not be changed.
  //
  // However, this function is also used to keep track (for statistical
  // purposes) of the vertical average of certain variables.  In that case,
  // the vertical average needs to be taken over the entire vertical domain
  // (level 1 to level gr%nz).
  //
  //
  // In both the thermodynamic-level computation and the momentum-level
  // computation, the numerator integral is divided by the denominator integral
  // in order to find the average value (over the vertical domain) of x.
  // References:
  // None
  //-----------------------------------------------------------------------
  // Input profile window: total_idx points, ordered from lowest to highest
  // altitude. rho_ds is dry static density [kg/m^3], field has arbitrary units,
  // and dz is level thickness [m]. The first elements are the caller's first
  // included level; no boundary levels are silently omitted here.
  static_assert(Rho::rank == 1 && Field::rank == 1 && Dz::rank == 1);
  const auto rho_ds_exec = to_execution_space(rho_ds);
  const auto field_exec = to_execution_space(field);
  const auto dz_exec = to_execution_space(dz);
  Real1D result("vertical_avg::result", 1);
  // A scalar result crosses to host once. The vertical sums stay ordered
  // inside one device work item rather than an unordered floating reduction.
  Kokkos::parallel_for(
      "vertical_avg", 1, KOKKOS_LAMBDA(const int) {
        real_t numer_integral = real_t(0.0);
        real_t denom_integral = real_t(0.0);
        for (int k = 0; k < total_idx; ++k) {
          numer_integral = numer_integral + rho_ds_exec(k) * dz_exec(k) * field_exec(k);
          denom_integral = denom_integral + rho_ds_exec(k) * dz_exec(k);
        }
        result(0) = numer_integral / denom_integral;
      });
  return to_host_space(result)(0);
}

// =============================================================================
template <class Rho, class Field, class Dz>
inline real_t vertical_integral(int total_idx, Rho rho_ds, // In
                                Field field, Dz dz         // In
) {
  // Description:
  // Computes the vertical integral. rho_ds, field, and dz must all be
  // of size total_idx and should all start at the same index.
  //
  // References:
  // None
  //-----------------------------------------------------------------------
  // Input profile window: total_idx points, ordered from lowest to highest
  // altitude. rho_ds is dry static density [kg/m^3], field has arbitrary units,
  // and dz is level thickness [m]. The first elements are the caller's first
  // included level; no boundary levels are silently omitted here.
  static_assert(Rho::rank == 1 && Field::rank == 1 && Dz::rank == 1);
  const auto rho_ds_exec = to_execution_space(rho_ds);
  const auto field_exec = to_execution_space(field);
  const auto dz_exec = to_execution_space(dz);
  Real1D result("vertical_integral::result", 1);
  // A scalar result crosses to host once. The vertical sums stay ordered
  // inside one device work item rather than an unordered floating reduction.
  Kokkos::parallel_for(
      "vertical_integral", 1, KOKKOS_LAMBDA(const int) {
        real_t vertical_integral = real_t(0.0);
        // Native SUM(field*rho_ds*dz) retains increasing profile order here.
        for (int k = 0; k < total_idx; ++k) {
          vertical_integral = vertical_integral + field_exec(k) * rho_ds_exec(k) * dz_exec(k);
        }
        result(0) = vertical_integral;
      });
  return to_host_space(result)(0);
}

// =============================================================================
inline void pvertinterp(int nzt, int ngrdcol, const grid& gr,         // In
                        Real2D p_mid, real_t p_out, Real2D input_var, // In
                        Real1D interp_var                             // Out
) {
  // Interpolate a column field to pressure p_out [Pa], extending the surface
  // or lid value outside its pressure range. Native bracket/search order and
  // k_upper=1 fallback are retained; no error is reported from unused l_error.

  // ------------------------ Input Variables ------------------------
  // p_mid: input level pressure levels
  // input_var: input  array
  // p_out: output pressure level
  // ------------------------ Output Variables ------------------------
  // interp_var: output array (interpolated)

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;

  // Initialize index array and logical flags
  // Native l_error is never consumed; keep it column-private below.
  // If we've fallen through the k=1,nz-1 loop, we cannot interpolate and
  // must extrapolate from the bottom or top data level for at least some
  // of the longitude points.
  Kokkos::parallel_for(
      "pvertinterp::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column + 1;
        int k{};
        int k_upper{};
        real_t dpu{};
        real_t dpl{};
        bool l_found{};
        [[maybe_unused]] bool l_error = false;
        if ((p_out >= p_mid(i - 1, gr_k_lb_zt - 1))) {
          interp_var(i - 1) = input_var(i - 1, gr_k_lb_zt - 1);
        } else if ((p_out <= p_mid(i - 1, gr_k_ub_zt - 1))) {
          interp_var(i - 1) = input_var(i - 1, gr_k_ub_zt - 1);
        } else {
          l_found = false;
          k_upper = 1;
          // Store level indices for interpolation.
          // If all indices for this level have been found,
          // do the interpolation
          for (k = gr_k_lb_zt; ((gr_grid_dir_indx) > 0 ? k <= ((gr_k_ub_zt - gr_grid_dir_indx))
                                                       : k >= ((gr_k_ub_zt - gr_grid_dir_indx)));
               k += gr_grid_dir_indx) {
            if (((p_mid(i - 1, k - 1) > p_out) &&
                 (p_out >= p_mid(i - 1, (k + gr_grid_dir_indx) - 1)))) {
              l_found = true;
              k_upper = k;
              break;
            }
          }
          if (!(l_found)) {
            l_error = true;
          }
          dpu = (p_mid(i - 1, k_upper - 1) - p_out);
          dpl = (p_out - p_mid(i - 1, (k_upper + gr_grid_dir_indx) - 1));
          interp_var(i - 1) = (((input_var(i - 1, k_upper - 1) * dpl) +
                                (input_var(i - 1, (k_upper + gr_grid_dir_indx) - 1) * dpu)) /
                               (dpl + dpu));
        }
      });
  return;
}

// =============================================================================
inline void calculate_thlp2_rad(int ngrdcol, int nzm, int nzt, const grid& gr,                 // In
                                Real2D rcm, Real2D thlprcp, Real2D radht, Real2D clubb_params, // In
                                Real2D thlp2_forcing // InOut
) {
  // Description:
  // Computes the contribution of radiative cooling to thlp2
  // References:
  // See clubb:ticket:632
  // ----------------------------------------------------------------------

  // Input Variables
  // ngrdcol: Number of grid columns
  // nzm: Number of momentum vertical levels              [-]
  // nzt: Number of thermodynamic vertical levels             [-]
  // rcm: Cloud water mixing ratio on thermodynamic grid      [kg/kg]
  // radht: SW + LW heating rate (on thermodynamic grid)        [K/s]
  // thlprcp: thl'rc'                                        [K kg/kg]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // Input/Output Variables
  // thlp2_forcing: <th_l'^2> forcing (momentum levels)            [K^2/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D rcm_zm("calculate_thlp2_rad::rcm_zm", ngrdcol,
                nzm); // Cloud water mixing ratio on momentum grid      [kg/kg]
  Real2D radht_zm("calculate_thlp2_rad::radht_zm", ngrdcol,
                  nzm); // SW + LW heating rate (on momentum grid)        [K/s]
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_nzm = gr.nzm;
  const auto gr_nzt = gr.nzt;

  Kokkos::deep_copy(rcm_zm, zt2zm_api(gr_nzm, gr_nzt, ngrdcol, gr, rcm));
  Kokkos::deep_copy(radht_zm, zt2zm_api(gr_nzm, gr_nzt, ngrdcol, gr, radht));
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calculate_thlp2_rad::cells_0", ngrdcol * level_count, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol + 1;
          const int k = cell / ngrdcol + 1;
          if ((rcm_zm(i - 1, k - 1) > rc_tol)) {
            thlp2_forcing(i - 1, k - 1) =
                (thlp2_forcing(i - 1, k - 1) +
                 ((((clubb_params(i - 1, ithlp2_rad_coef) * two) * radht_zm(i - 1, k - 1)) /
                   rcm_zm(i - 1, k - 1)) *
                  thlprcp(i - 1, k - 1)));
          }
        });
  }
  return;
}

// =============================================================================
// Native generic dispatch. These overloads contain no physics or extra options.
inline void calc_xpwp(const grid& gr, Real1D Km_zm, Real1D xm, Real1D xpwp) {
  // Select the same-ranked Native procedure.
  calc_xpwp_1D(gr, Km_zm, xm, xpwp);
}

inline void calc_xpwp(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D Km_zm, Real2D xm,
                      Real2D xpwp) {
  // Select the same-ranked Native procedure.
  calc_xpwp_2D(nzm, nzt, ngrdcol, gr, Km_zm, xm, xpwp);
}

inline Real2D smooth_min(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_min_sclr_idx(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

inline Real2D smooth_min(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_min_array_scalar(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

inline Real2D smooth_min(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_min_arrays(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

KOKKOS_INLINE_FUNCTION real_t smooth_min(real_t input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_min_scalars(input_var1, input_var2, smth_coef);
}

inline Real2D smooth_max(int nz, int ngrdcol, real_t input_var1, Real2D input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_max_sclr_idx(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

inline Real2D smooth_max(int nz, int ngrdcol, Real2D input_var1, real_t input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_max_array_scalar(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

inline Real2D smooth_max(int nz, int ngrdcol, Real2D input_var1, Real2D input_var2,
                         real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_max_arrays(nz, ngrdcol, input_var1, input_var2, smth_coef);
}

inline Real1D smooth_max(int ngrdcol, Real1D input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_max_array_1D_scalar(ngrdcol, input_var1, input_var2, smth_coef);
}

KOKKOS_INLINE_FUNCTION real_t smooth_max(real_t input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same-ranked Native procedure.
  return smooth_max_scalars(input_var1, input_var2, smth_coef);
}

// =============================================================================
// Host/strided adapters stage only at the existing public array API boundaries.
// Conversions use array_types; core execution-space calls select the exact
// overloads above and keep their allocation throughout the call chain.

template <class Input>
inline std::enable_if_t<Kokkos::is_view<Input>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space>>
smooth_min(int nz, int ngrdcol, real_t input_var1, Input input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input::rank == 2);
  const auto result =
      smooth_min_sclr_idx(nz, ngrdcol, input_var1, to_execution_space(input_var2), smth_coef);
  if constexpr (std::is_same_v<typename Input::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space> output(
        "smooth_min::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input>
inline std::enable_if_t<Kokkos::is_view<Input>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space>>
smooth_min(int nz, int ngrdcol, Input input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input::rank == 2);
  const auto result =
      smooth_min_array_scalar(nz, ngrdcol, to_execution_space(input_var1), input_var2, smth_coef);
  if constexpr (std::is_same_v<typename Input::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space> output(
        "smooth_min::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input1, class Input2>
inline std::enable_if_t<Kokkos::is_view<Input1>::value && Kokkos::is_view<Input2>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input1::memory_space>>
smooth_min(int nz, int ngrdcol, Input1 input_var1, Input2 input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input1::rank == 2 && Input2::rank == 2);
  const auto result = smooth_min_arrays(nz, ngrdcol, to_execution_space(input_var1),
                                        to_execution_space(input_var2), smth_coef);
  if constexpr (std::is_same_v<typename Input1::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input1::memory_space> output(
        "smooth_min::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input>
inline std::enable_if_t<Kokkos::is_view<Input>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space>>
smooth_max(int nz, int ngrdcol, real_t input_var1, Input input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input::rank == 2);
  const auto result =
      smooth_max_sclr_idx(nz, ngrdcol, input_var1, to_execution_space(input_var2), smth_coef);
  if constexpr (std::is_same_v<typename Input::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space> output(
        "smooth_max::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input>
inline std::enable_if_t<Kokkos::is_view<Input>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space>>
smooth_max(int nz, int ngrdcol, Input input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input::rank == 2);
  const auto result =
      smooth_max_array_scalar(nz, ngrdcol, to_execution_space(input_var1), input_var2, smth_coef);
  if constexpr (std::is_same_v<typename Input::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input::memory_space> output(
        "smooth_max::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input1, class Input2>
inline std::enable_if_t<Kokkos::is_view<Input1>::value && Kokkos::is_view<Input2>::value,
                        Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input1::memory_space>>
smooth_max(int nz, int ngrdcol, Input1 input_var1, Input2 input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input1::rank == 2 && Input2::rank == 2);
  const auto result = smooth_max_arrays(nz, ngrdcol, to_execution_space(input_var1),
                                        to_execution_space(input_var2), smth_coef);
  if constexpr (std::is_same_v<typename Input1::memory_space, typename Real2D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t**, Kokkos::LayoutLeft, typename Input1::memory_space> output(
        "smooth_max::output", ngrdcol, nz);
    from_execution_space(output, result);
    return output;
  }
}

template <class Input>
inline std::enable_if_t<Kokkos::is_view<Input>::value,
                        Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Input::memory_space>>
smooth_max(int ngrdcol, Input input_var1, real_t input_var2, real_t smth_coef) {
  // Select the same Native rank and return owned storage in the input's space.
  static_assert(Input::rank == 1);
  const auto result =
      smooth_max_array_1D_scalar(ngrdcol, to_execution_space(input_var1), input_var2, smth_coef);
  if constexpr (std::is_same_v<typename Input::memory_space, typename Real1D::memory_space>) {
    return result;
  } else {
    Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Input::memory_space> output(
        "smooth_max::output", ngrdcol);
    from_execution_space(output, result);
    return output;
  }
}

template <class Lhs>
inline void
set_boundary_conditions_lhs(int diag_index, int low_bound, int high_bound, // In
                            Lhs lhs,                                       // InOut
                            std::optional<int> diag_index2 = std::nullopt, // Optional In
                            std::optional<int> low_bound2 = std::nullopt,  // Optional In
                            std::optional<int> high_bound2 = std::nullopt  // Optional In
) {
  // Preserve caller-owned host/strided LAPACK storage and Native indices.
  static_assert(Lhs::rank == 2);
  auto lhs_exec = to_execution_space(lhs);
  set_boundary_conditions_lhs(diag_index, low_bound, high_bound, lhs_exec, diag_index2, low_bound2,
                              high_bound2);
  from_execution_space(lhs, lhs_exec);
}

template <class Rhs>
inline void
set_boundary_conditions_rhs(real_t low_value, int low_bound, real_t high_value,
                            int high_bound,                                   // In
                            Rhs rhs,                                          // InOut
                            std::optional<real_t> low_value2 = std::nullopt,  // Optional In
                            std::optional<int> low_bound2 = std::nullopt,     // Optional In
                            std::optional<real_t> high_value2 = std::nullopt, // Optional In
                            std::optional<int> high_bound2 = std::nullopt     // Optional In
) {
  // Preserve caller-owned host/strided LAPACK storage and overwrite order.
  static_assert(Rhs::rank == 1);
  auto rhs_exec = to_execution_space(rhs);
  set_boundary_conditions_rhs(low_value, low_bound, high_value, high_bound, rhs_exec, low_value2,
                              low_bound2, high_value2, high_bound2);
  from_execution_space(rhs, rhs_exec);
}

template <class Km, class Xm, class Xpwp>
inline void calc_xpwp(const grid& gr, Km Km_zm, Xm xm, Xpwp xpwp) {
  // Stage the Native 1D generic; only interior output levels are overwritten.
  static_assert(Km::rank == 1 && Xm::rank == 1 && Xpwp::rank == 1);
  const auto Km_exec = to_execution_space(Km_zm);
  const auto xm_exec = to_execution_space(xm);
  auto xpwp_exec = to_execution_space(xpwp);
  calc_xpwp_1D(gr, Km_exec, xm_exec, xpwp_exec);
  from_execution_space(xpwp, xpwp_exec);
}

template <class Km, class Xm, class Xpwp>
inline void calc_xpwp(int nzm, int nzt, int ngrdcol, const grid& gr, Km Km_zm, Xm xm, Xpwp xpwp) {
  // Stage the Native 2D generic, retaining caller-owned boundary storage.
  static_assert(Km::rank == 2 && Xm::rank == 2 && Xpwp::rank == 2);
  const auto Km_exec = to_execution_space(Km_zm);
  const auto xm_exec = to_execution_space(xm);
  auto xpwp_exec = to_execution_space(xpwp);
  calc_xpwp_2D(nzm, nzt, ngrdcol, gr, Km_exec, xm_exec, xpwp_exec);
  from_execution_space(xpwp, xpwp_exec);
}

} // namespace clubb_kokkos
