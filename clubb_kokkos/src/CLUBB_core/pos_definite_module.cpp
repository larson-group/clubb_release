// Limit scalar fluxes while retaining their conservative field update.
//
// pos_definite_adj changes field_np1 and flux_np1 in place and fills field_pd
// and flux_pd budgets. Field arrays are (column,thermodynamic level), flux
// arrays are (column,momentum level), all LayoutLeft. The physical-direction
// limiter loop and storage-endpoint fluxes follow Native grid metadata. Two
// separate column ANY checks select budgets before/after the field update.
// All scratch is owned per call.
#include "scratch_workspace.hpp"
#include "error_code.hpp"
#include "pos_definite_module.hpp"
#include "constants_clubb.hpp"
#include <iostream>

namespace clubb_kokkos {

//-----------------------------------------------------------------------
// dt: Timestep [s].
// field_np1: Field at n+1 (e.g. rtm [kg/kg]); updated in place.
// flux_np1: Flux applied to the field; replaced with limited flux.
// field_n: Field at n, prior to n+1, with the same units as field_np1.
// field_pd: Budget of the change in the mean field due to the scheme [field units/s].
// flux_pd: Budget of the flux change due to the scheme [flux units/s].
void pos_definite_adj(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                      real_t dt,                         // In
                      Real2D field_np1, Real2D flux_np1, // InOut
                      Real2D field_n,                    // In
                      Real2D field_pd, Real2D flux_pd)   // Out
{
  // Description:
  //   Applies a  flux conservative positive definite scheme to a variable

  //   There are two possible grids:
  //   (1) flux on zm  field on zt
  //   then
  //   flux_zt(k) = ( flux_zm(k+1) + flux_zm(k) ) / 2

  //         CLUBB grid                  Smolarkiewicz grid
  //   m +-- flux  zm(k+1) --+               flux        k + 1/2
  //   t +-- field zt(k)   --+               field, fout k
  //   m +-- flux  zm(k)   --+               flux        k - 1/2
  //   t +-- field zt(k-1) --+

  //   (2) flux on zt field on zm
  //   then
  //   flux_zm(k) = ( flux_zt(k) + flux_zt(k-1) ) / 2

  //         CLUBB grid                  Smolarkiewicz grid
  //   m +-- field  (k+1)  --+
  //   t +-- flux   (k)    --+               flux        k + 1/2
  //   m +-- field  (k)    --+               field, fout k
  //   t +-- flux   (k-1)  --+               flux        k - 1/2

  // References:
  //   ``A Positive Definite Advection Scheme Obtained by
  //     Nonlinear Renormalization of the Advective Fluxes'' Smolarkiewicz (1989)
  //     Monthly Weather Review, Vol. 117, pp. 2626--2632
  //-----------------------------------------------------------------------

  // Native local arrays are independent, owned per-call allocations.
  // [F_i+1/2]^+ and [F_i+1/2]^- in Smolarkiewicz [flux units].
  auto workspace_scope = workspace.scope();
  Real2D flux_plus     = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "pos_definite_adj::flux_plus", ngrdcol, nzm);
  Real2D flux_minus    = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "pos_definite_adj::flux_minus", ngrdcol, nzm);
  // Limited flux at n+1 [flux units].
  Real2D flux_lim = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                 "pos_definite_adj::flux_lim", ngrdcol, nzm);
  // Eqn A4: F_i^OUT, total outward scalar-level flux [flux units].
  Real2D fout = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "pos_definite_adj::fout",
                                             ngrdcol, nzt);
  // Original field at n+1, retained for the budget [field units].
  Real2D field_nonlim = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "pos_definite_adj::field_nonlim", ngrdcol, nzt);
  // Grid-spacing/timestep conversion factor [m/s].
  Real2D dz_over_dt = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                   "pos_definite_adj::dz_over_dt", ngrdcol, nzt);

  // Capture grid components; boundary metadata uses zero-based slots.
  const auto invrs_dzt    = gr.invrs_dzt;
  const real_t grid_dir   = gr.grid_dir;
  const int grid_dir_indx = gr.grid_dir_indx;
  using Policy            = Kokkos::RangePolicy<typename Real2D::execution_space>;

  //-----------------------------------------------------------------------
  // TODO: Replace this host-only debug query when its Native owner is ported.
  if (clubb_at_least_debug_level_api(1)) {
    std::cout << "Correcting flux\n";
  }

  Kokkos::parallel_for(
      "pos_definite_adj::flux_split", Policy(0, ngrdcol * nzm), KOKKOS_LAMBDA(int cell) {
        const int k = cell / ngrdcol;
        const int i = cell % ngrdcol;
        // Def. of F+ and F- from eqn 2 Smolarkowicz
        // Materialize same-kind constant values for device-callable min/max.
        flux_plus(i, k)  = Kokkos::max(real_t(zero_threshold), flux_np1(i, k));
        flux_minus(i, k) = -Kokkos::min(real_t(zero_threshold), flux_np1(i, k));
      });

  Kokkos::parallel_for(
      "pos_definite_adj::outward_flux", Policy(0, ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int k = cell / ngrdcol;
        const int i = cell % ngrdcol;
        // These independent scalar-level fields share the same cell traversal.
        dz_over_dt(i, k) = (real_t(1.) / (grid_dir * invrs_dzt(i, k))) / dt;
        int kphalf, kmhalf;
        // If the scalar variable is on the kth t-level, then
        // Smolarkowicz's k+1/2 flux level is the k+1 m-level in CLUBB.
        if (grid_dir_indx > 0) {
          // Ascending grid
          kphalf = Kokkos::min(k + 1, nzm - 1); // k+1/2 flux level
          kmhalf = Kokkos::max(k, 0);           // k-1/2 flux level
        } else {                                // grid_dir_indx < 0
          // Descending grid
          kphalf = Kokkos::max(k, 0);           // k+1/2 flux level
          kmhalf = Kokkos::min(k + 1, nzm - 1); // k-1/2 flux level
        }

        // Eqn A4 from Smolarkowicz
        // We place a limiter of eps to prevent a divide by zero, and
        //   after this calculation fout is on the scalar level, and
        //   fout is the total outward flux for the scalar level k.
        fout(i, k) = Kokkos::max(flux_plus(i, kphalf) + flux_minus(i, kmhalf), real_t(eps));
      });

  // Flatten the original physical-direction DO by ordinal, retaining its exact
  // Native bounds/step. Each iteration writes a distinct interior flux level;
  // this loop has no recurrence. Storage endpoint copies remain below.
  const int first_k     = gr.k_lb_zm;
  const int last_k      = gr.k_ub_zm - 2 * grid_dir_indx;
  const int level_count = Kokkos::max(0, (last_k - first_k) / grid_dir_indx + 1);
  Kokkos::parallel_for(
      "pos_definite_adj::limit", Policy(0, ngrdcol * level_count), KOKKOS_LAMBDA(int cell) {
        const int k = first_k + (cell / ngrdcol) * grid_dir_indx;
        const int i = cell % ngrdcol;
        int kphalf, kp1, k_zt;
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // FIXME:
        // We haven't tested this for negative values at the nz level
        // -dschanen 13 June 2008
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (grid_dir_indx > 0) {
          // Ascending grid
          kphalf = Kokkos::min(k + 1, nzm - 1); // k+1/2 flux level
          kp1    = Kokkos::min(k + 1, nzt - 1); // k+1 scalar level
          k_zt   = Kokkos::min(k, nzt - 1);     // k scalar level
        } else {                                // grid_dir_indx
          kphalf = Kokkos::max(k - 1, 0);       // k+1/2 flux level
          kp1    = Kokkos::max(k - 2, 0);       // k+1 scalar level
          k_zt   = Kokkos::max(k - 1, 0);       // k scalar level
        }

        // Eqn 10 from Smolarkowicz (1989)
        flux_lim(i, kphalf) = Kokkos::max(
            Kokkos::min(flux_np1(i, kphalf), (flux_plus(i, kphalf) / fout(i, k_zt)) *
                                                 field_n(i, k_zt) * dz_over_dt(i, k_zt)),
            -((flux_minus(i, kphalf) / fout(i, kp1)) * field_n(i, kp1) * dz_over_dt(i, k_zt)));
      });

  // Boundary conditions
  Kokkos::parallel_for(
      "pos_definite_adj::boundaries", Policy(0, ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int i          = cell;
        flux_lim(i, 0)       = flux_np1(i, 0);
        flux_lim(i, nzm - 1) = flux_np1(i, nzm - 1);
      });

  Kokkos::parallel_for(
      "pos_definite_adj::flux_budget", Policy(0, ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int i = cell;
        // Only set flux_pd for a column if there is a below zero value in that column.
        // Keep this Native ANY before the field update; no host staging is needed.
        bool l_negative = false;
        for (int k = 0; k < nzt; ++k) {
          l_negative = l_negative || field_np1(i, k) < zero;
        }
        if (l_negative) {
          for (int k = 0; k < nzm; ++k) {
            flux_pd(i, k) = (flux_lim(i, k) - flux_np1(i, k)) / dt;
          }
        } else {
          for (int k = 0; k < nzm; ++k) {
            flux_pd(i, k) = zero;
          }
        }
      });

  // Native assignment copies values; a View alias would lose the original field.
  Kokkos::deep_copy(field_nonlim, field_np1);

  // Apply change to field at n+1
  // Materialize the Native array-expression argument at the existing ddzm API.
  Real2D flux_correction = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "pos_definite_adj::flux_correction", ngrdcol, nzm);
  Kokkos::parallel_for(
      "pos_definite_adj::flux_correction", Policy(0, ngrdcol * nzm), KOKKOS_LAMBDA(int cell) {
        const int k           = cell / ngrdcol;
        const int i           = cell % ngrdcol;
        flux_correction(i, k) = flux_lim(i, k) - flux_np1(i, k);
      });
  Real2D ddzm_flux_correction = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "pos_definite_adj::ddzm_flux_correction", ngrdcol, nzt);
  ddzm(nzm, nzt, ngrdcol, gr, flux_correction, ddzm_flux_correction);
  Kokkos::parallel_for(
      "pos_definite_adj::field_update", Policy(0, ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int k     = cell / ngrdcol;
        const int i     = cell % ngrdcol;
        field_np1(i, k) = -dt * ddzm_flux_correction(i, k) + field_np1(i, k);
      });

  // Determine the total time tendency in field due to this calculation
  // (for diagnostic purposes)
  Kokkos::parallel_for(
      "pos_definite_adj::field_budget", Policy(0, ngrdcol), KOKKOS_LAMBDA(int cell) {
        const int i = cell;
        // Only set field_pd for a column if there is a below zero value in that column.
        // This second Native ANY examines the updated field, not the earlier mask.
        bool l_negative = false;
        for (int k = 0; k < nzt; ++k) {
          l_negative = l_negative || field_np1(i, k) < zero;
        }
        if (l_negative) {
          for (int k = 0; k < nzt; ++k) {
            field_pd(i, k) = (field_np1(i, k) - field_nonlim(i, k)) / dt;
          }
        } else {
          for (int k = 0; k < nzt; ++k) {
            field_pd(i, k) = zero;
          }
        }
      });

  // Replace the non-limited flux with the limited flux
  Kokkos::deep_copy(flux_np1, flux_lim);
}

//-----------------------------------------------------------------------
} // namespace clubb_kokkos
