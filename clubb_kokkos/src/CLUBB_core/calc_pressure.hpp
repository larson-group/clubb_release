// Compute hydrostatic pressure and linearized mean virtual potential temperature.
//
// init_pressure fills pressure [Pa] and Exner[-] on both staggered grids from
// theta_v[K] and surface pressure[Pa]. calculate_thvm fills theta_v[K] from
// liquid-water potential temperature, total/cloud water and a reference theta_v.
// Both preserve Native argument order and column-first fields. The hydrostatic
// recurrence follows increasing storage levels within each independent column;
// its surface is storage level 0 (Native level 1). Scratch uses the caller workspace.
// Host/strided init callers use the existing array owner; calculate_thvm runs
// directly in its output View's execution space.
#pragma once
#include "array_types.hpp"
#include "scratch_workspace.hpp"
#include "constants_clubb.hpp"
#include "grid_class.hpp"
#include <limits>

namespace clubb_kokkos {

//=============================================================================
inline void init_pressure(ScratchWorkspace& workspace, int ngrdcol, const grid& gr, Real2D thvm,
                          Real1D p_sfc, // In
                          Real2D p_in_Pa, Real2D exner, Real2D p_in_Pa_zm,
                          Real2D exner_zm) // Out
{
  // Description:
  // Calculates the initial pressure according to the hydrostatic
  // approximation.  Combining the moist equation of state and the hydrostatic
  // approximation, the change of pressure with respect to height can be
  // calculated based on theta_v, such that:
  //
  // dp/dz = - p * grav / ( Rd * theta_v * exner );
  //
  // where exner = ( p / p0 )^(Rd/Cp);
  //
  // and where p0 is a reference pressure of 100000 Pa.
  //
  // The integral equation is set up to integrate over p on the left-hand side
  // and integrate over z on the right-hand side.  The equation is:
  //
  // INT(p1:p2) p^(Rd/Cp-1) dp
  // = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.
  //
  // The value of mean theta_v (thvm) is calculated at each thermodynamic grid
  // level, and linear interpolation is used in the integral equation for all
  // altitude in-between successive thermodynamic levels, such that:
  //
  // thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
  //
  // The integrals are solved, and the results for pressure can be rewritten
  // in terms of exner, such that:
  //
  // exner2 - exner1
  //   | - ( grav / Cp )
  //   |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
  // = | where thvm2 /= thvm1;
  //   |
  //   | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
  //
  // The value of pressure (exner) can be calculated using the above equation
  // at all vertical levels once the value of pressure is known at one level.
  // Since the surface pressure is known at the initial time, that allows
  // pressure to be calculated for the rest of the vertical profile.

  // References:
  //-----------------------------------------------------------------------

  // thvm: Mean theta_v on thermodynamic levels[K]. p_sfc: surface pressure[Pa].
  // p_in_Pa/exner: (column,thermodynamic level)[Pa]/[-].
  // p_in_Pa_zm/exner_zm: (column,momentum level)[Pa]/[-].
  constexpr real_t g_ov_Cp  = grav / Cp; // g/Cp[K/m], Native constant parameter.
  const real_t thvm_epsilon = std::numeric_limits<real_t>::epsilon();
  // Capture component Views/scalars; owning grid metadata stays on the host.
  const auto gr_zt = gr.zt;
  const auto gr_zm = gr.zm;
  const int nzt    = gr.nzt;
  const int nzm    = gr.nzm;
  using Policy     = Kokkos::RangePolicy<typename Real2D::execution_space>;

  // The pressure (and Exner) at the first momentum storage level is the surface
  // pressure (and Exner based on the surface pressure), as explicitly in Native.
  Kokkos::parallel_for(
      "init_pressure::surface", Policy(0, ngrdcol), KOKKOS_LAMBDA(int column) {
        const int i      = column;
        p_in_Pa_zm(i, 0) = p_sfc(i);
        exner_zm(i, 0)   = Kokkos::pow(p_sfc(i) / p0, kappa);
      });

  // Interpolate theta_v into this call's scoped momentum-level scratch.
  auto workspace_scope = workspace.scope();
  const auto thvm_zm   = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "init_pressure::thvm_zm", ngrdcol, nzm);
  zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, thvm, thvm_zm, zero_threshold);

  // Calculate Exner at all other thermodynamic and momentum grid levels.
  // exner2
  // = exner1
  //     | ( grav / Cp )
  //     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
  //   - | where thvm2 /= thvm1;
  //     |
  //     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).

  // Calculate Exner at thermodynamic storage level 0 (Native 1). Do not reinterpret the
  // source's fixed storage endpoint using physical grid boundary metadata.
  Kokkos::parallel_for(
      "init_pressure::first_exner", Policy(0, ngrdcol), KOKKOS_LAMBDA(int column) {
        const int i = column;
        if (Kokkos::abs(thvm(i, 0) - thvm_zm(i, 0)) > thvm_epsilon * thvm(i, 0)) {
          exner(i, 0) = exner_zm(i, 0) - g_ov_Cp * (gr_zt(i, 0) - gr_zm(i, 0)) /
                                             (thvm(i, 0) - thvm_zm(i, 0)) *
                                             Kokkos::log(thvm(i, 0) / thvm_zm(i, 0));
        } else { // thvm(1) = thvm_zm(1).
          exner(i, 0) = exner_zm(i, 0) - g_ov_Cp * (gr_zt(i, 0) - gr_zm(i, 0)) / thvm(i, 0);
        }
      });

  // Calculate pressure on thermodynamic levels.
  Kokkos::parallel_for(
      "init_pressure::first_pressure", Policy(0, ngrdcol), KOKKOS_LAMBDA(int column) {
        const int i   = column;
        p_in_Pa(i, 0) = p0 * Kokkos::pow(exner(i, 0), one / kappa);
      });

  // Loop over all other thermodynamic vertical grid levels. Each independent
  // column keeps the Native k=2:nzt recurrence and per-level calculation order.
  Kokkos::parallel_for(
      "init_pressure::thermodynamic", Policy(0, ngrdcol), KOKKOS_LAMBDA(int column) {
        const int i = column;
        for (int k = 1; k < nzt; ++k) {
          // Calculate Exner on thermodynamic levels.
          if (Kokkos::abs(thvm(i, k) - thvm(i, k - 1)) > thvm_epsilon * thvm(i, k)) {
            exner(i, k) = exner(i, k - 1) - g_ov_Cp * (gr_zt(i, k) - gr_zt(i, k - 1)) /
                                                (thvm(i, k) - thvm(i, k - 1)) *
                                                Kokkos::log(thvm(i, k) / thvm(i, k - 1));
          } else { // Adjacent thermodynamic theta_v values agree at the Native threshold.
            exner(i, k) = exner(i, k - 1) - g_ov_Cp * (gr_zt(i, k) - gr_zt(i, k - 1)) / thvm(i, k);
          }
          // Calculate pressure on thermodynamic levels.
          p_in_Pa(i, k) = p0 * Kokkos::pow(exner(i, k), one / kappa);
        }
      });

  // Loop over all momentum grid levels. Each remaining momentum cell depends
  // on thermo(k-1), not on a preceding momentum update, so these are independent.
  Kokkos::parallel_for(
      "init_pressure::momentum", Policy(0, ngrdcol * (nzm - 1)), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol + 1;
        // Calculate Exner on momentum levels. Keep epsilon*thvm(k-1), as in Native.
        if (Kokkos::abs(thvm(i, k - 1) - thvm_zm(i, k)) > thvm_epsilon * thvm(i, k - 1)) {
          exner_zm(i, k) = exner(i, k - 1) - g_ov_Cp * (gr_zm(i, k) - gr_zt(i, k - 1)) /
                                                 (thvm_zm(i, k) - thvm(i, k - 1)) *
                                                 Kokkos::log(thvm_zm(i, k) / thvm(i, k - 1));
        } else { // Adjacent thermo/momentum theta_v values agree at the Native threshold.
          exner_zm(i, k) =
              exner(i, k - 1) - g_ov_Cp * (gr_zm(i, k) - gr_zt(i, k - 1)) / thvm_zm(i, k);
        }
        // Calculate pressure on momentum levels.
        p_in_Pa_zm(i, k) = p0 * Kokkos::pow(exner_zm(i, k), one / kappa);
      });
}

//=============================================================================
template <class Thlm, class Rtm, class Rcm, class Exner, class ThvDsZt, class Thvm>
void calculate_thvm(int nzt, int ngrdcol,                                        // In
                    Thlm thlm, Rtm rtm, Rcm rcm, Exner exner, ThvDsZt thv_ds_zt, // In
                    Thvm thvm)                                                   // Out
{
  // Description:
  // Calculates mean theta_v based on a linearized approximation to the theta_v
  // equation.  This equation also includes liquid water loading.

  // References:
  //-----------------------------------------------------------------------

  // thlm: Mean theta_l[K]; rtm/rcm: total/cloud water[kg/kg]; exner[-];
  // thv_ds_zt: reference theta_v[K]. All fields are(column,thermodynamic level).
  static_assert(Thlm::rank == 2 && Rtm::rank == 2 && Rcm::rank == 2 && Exner::rank == 2 &&
                ThvDsZt::rank == 2 && Thvm::rank == 2);
  // Keep host/strided driver work on the host and core work on its execution
  // space. Input Views must be accessible from the caller output's space.
  using Policy = Kokkos::RangePolicy<typename Thvm::execution_space>;
  // Calculate mean theta_v.
  Kokkos::parallel_for(
      "calculate_thvm", Policy(0, ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        thvm(i, k)  = thlm(i, k) + ep1 * thv_ds_zt(i, k) * rtm(i, k) +
                      (Lv / (Cp * exner(i, k)) - ep2 * thv_ds_zt(i, k)) * rcm(i, k);
      });
}

// Preserve host/strided initialization callers at the existing shared boundary.
template <class Thvm, class Surface, class Pressure, class Exner, class MomentumPressure,
          class MomentumExner>
void init_pressure(ScratchWorkspace& workspace, int ngrdcol, const grid& gr, Thvm thvm,
                   Surface p_sfc, // In
                   Pressure p_in_Pa, Exner exner, MomentumPressure p_in_Pa_zm,
                   MomentumExner exner_zm) // Out
{
  // Grid components already live in the core execution space; stage once here.
  auto workspace_scope = workspace.scope();
  auto thvm_exec       = to_execution_space(workspace, thvm);
  auto p_sfc_exec      = to_execution_space(workspace, p_sfc);
  auto p_exec          = to_execution_space(workspace, p_in_Pa, false);
  auto exner_exec      = to_execution_space(workspace, exner, false);
  auto p_zm_exec       = to_execution_space(workspace, p_in_Pa_zm, false);
  auto exner_zm_exec   = to_execution_space(workspace, exner_zm, false);
  init_pressure(workspace, ngrdcol, gr, thvm_exec, p_sfc_exec, // In
                p_exec, exner_exec, p_zm_exec, exner_zm_exec); // Out
  from_execution_space(workspace, p_in_Pa, p_exec);
  from_execution_space(workspace, exner, exner_exec);
  from_execution_space(workspace, p_in_Pa_zm, p_zm_exec);
  from_execution_space(workspace, exner_zm, exner_zm_exec);
}
} // namespace clubb_kokkos
