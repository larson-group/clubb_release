// Simplified radiation; see the matching public header for scheme ownership.
// Historical COAMPS/CLUBB grid context retained from the Native module.
// The simplified LW radiation code used by both various idealized cases is now
// consolidated into the simple_rad subroutine below.  The comments here are
// in regard to the differences between the COAMPS-LES simulations done in the
// past and how the idealized radiation is handled here in CLUBB.
// -dschanen 29 Sep 2010

// SPECIAL METHOD USED TO CALCULATE RADIATION
// Grid descriptions by Adam Smith, 27 June 2006
//
// In order to verify our CLUBB simulations are working properly, we
// have first developed a series of 3D simulations using the COAMPS-LES
// model.  This large-eddy simulation (LES) simulation uses specific
// methods to calculate radiation, subsidence, and other microphysical
// processes.  To make the two models simluate clouds as closely as
// possible, we use the same radiation scheme in both models.
//
// In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement
// all radiation code.  This allows the subroutine to be duplicated
// exactly in many different models.  However, the subroutine uses the
// COAMPS vertical grid.  Therefore, for CLUBB to implement this code
// correctly, we must modify some of our variable profiles before
// calling the radiation subroutine.
//
// The following diagram describes the differences in model grids
// (see "ADDITIONAL NOTE" for important grid information):
//
//       COAMPS-LES                                   CLUBB
//
//  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)
//
//  k= 1   (t) ----------                   ---------- k=kk        (t)
//
//  k= 2   (m) ----------                   ---------- k=kk        (m)
//
//  k= 2   (t) ----------                   ---------- k=kk-1      (t)
//
//            .                  .                     .
//            .                  .                     .
//            .                  .                     .
//
//  k=kk-1 (m) ----------  m = momentum     ---------- k=3         (m)
//                                  levels
//  k=kk-1 (t) ----------  t = thermo       ---------- k=2         (t)
//                                  levels
//  k=kk   (m) ----------                   ---------- k=2         (m)
//
//  k=kk   (t) ----------  kk = number of   ---------- k=1         (t)
//                              vertical
//  k=kk+1 (m) ----------       heights     ---------- k=1         (m)
//
// //////////////////////// MODEL SURFACE /////////////////////////////
//
// ADDITIONAL NOTE: In order to reconcile the COAMPS and CLUBB grids,
//                  the uppermost level of the COAMPS grid is ignored,
//                  eliminating the need to add level kk+1 to the CLUBB
//                  grid.  For the purposes of this code, the COAMPS
//                  grid levels are re-indexed to start from level 1
//                  at the uppermost useful COAMPS grid level (which
//                  was previously referred to as level 2 in the above
//                  diagram and is at the same altitude as CLUBB
//                  level kk).  Likewise, the lowermost COAMPS grid
//                  level is now indexed at kk (rather than kk+1 in the
//                  above diagram).  Brian Griffin; May 10, 2008.
//
// Also, the COAMPS grid indices are numbered from the top of the model
// downward, while the CLUBB grid indices are numbered from the bottom
// up.  Therefore, since we are using a COAMPS radiation scheme, we
// flip moisture and temperature profiles that are passed into the
// rad_lwsw subroutine.  The rad scheme will produce results in using
// the COAMPS grid scheme, so all radiation output will be flipped
// back to the CLUBB grid before being applied to the model.
//
// Finally, since the COAMPS scheme does not have a gridpoint below
// model surface, we add that point to all radiative output files once
// they are converted back to CLUBB setup.  This allows all averages and
// calculations to be done correctly.
//
//
// Computation of radiative fluxes on staggered grid
// Comments by Michael Falk, 16 February 2005.
//
// Frad (and its components Frad_LW and Frad_SW) should be computed on
// w points, not on mass points, which is apparent from its formulation
// and from its location in stats_sw instead of stats_sm.  The grid
// looks like this:
//
//
// -----Frad----------------------------------    k = 1  (w level)
//     /    \            |-dwm
// -LWP------radht----------------------------    k = 1  (mass level)
//     \    /            |
// -----Frad----------------------------------    k = 2  (w level)
//     /    \
// -LWP------radht----------------------------    k = 2  (mass level)
//     \    /
// -----Frad----------------------------------    k = 3  (w level)
//     /    \
// -LWP------radht----------------------------    k = 3  (mass level)
//
// If you consider Frad to take place on mass levels, then computing
// LWP is a forward difference and is only first-order accurate, while
// if Frad computed in between LWP levels, it is a centered difference
// which is second-order accurate.
//
// The coding implementation requires that Frad depend on LWP(k) and
// LWP(k-1) since the w level for a given k is at a higher altitude
// than the mass level.  radht, back on mass levels, depends on Frad(k)
// and Frad(k+1).
//
// ADDITIONAL NOTE: For clarification of terminology, a w level on the
//                  COAMPS grid is equivalent to a momentum (m) level
//                  on the CLUBB grid, and a mass level on the COAMPS
//                  grid is equivalent to a thermodynamic (t) level on
//                  the CLUBB grid.  Brian Griffin; May 10, 2008.
//
// Additionally, these computations assume that the distance between
// mass levels (dsigma) is constant, and that the w levels (spaced by
// dsigmw) always fall exactly halfway in between the mass levels.  If
// this is not the case, consider dwm to be the distance between a w
// level and the mass level below it, and dmw to be the distance
// between a mass level and the w level below it.  Then, the
// formulation for Frad_LW, for instance, would use a weighted average:
//
// (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
// which, for dwm always == dmw, reduces to
// (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
// which is identical to the current formulation.
// ((lwp(k)+lwp(k-1))/2)
//
// ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
//                  compatible with the use of a stretched
//                  (or unevenly-spaced) grid, as well as with the use
//                  of an evenly-spaced grid.  Interpolation functions
//                  are used to compute any weighted averages, rather
//                  than using general numbers such as (1/2), which is
//                  compatible only with an evenly-spaced grid.
//                  Brian Griffin; May 10, 2008.
//
//
//-----------------------------------------------------------------------

#include "scratch_workspace.hpp"
#include "CLUBB_core/file_functions.hpp"
#include "simple_rad_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include "Radiation/radiation_io.hpp"

namespace clubb_kokkos {
namespace {
constexpr int lba_ntimes = 36;
constexpr int lba_nzrad  = 33;
Real2D liq_water_path(ScratchWorkspace& workspace,               // InOut
                      int ngrdcol, int nzm, int nzt,             // In
                      Real2D rho, Real2D rcm, Real2D invrs_dzt); // In
} // namespace

void simple_rad(ScratchWorkspace& workspace,                            // InOut
                const grid& gr, int ngrdcol, Real2D rho, Real2D rho_zm, // In
                Real2D rtm, Real2D rcm, Real2D exner,                   // In
                stats_type& stats,                                      // InOut
                Real2D Frad_LW, Real2D radht_LW,                        // Out
                const parameters_radiation_state& parameters_radiation  // Native module In
) {
  // Compute simplified longwave flux and heating. Liquid water path and flux
  // lie on momentum levels; heating is their thermodynamic-level divergence.
  // F0 and F1 describe cloud-top/base heating, kappa absorption [m^2/kg].

  auto workspace_scope    = workspace.scope();
  constexpr real_t ls_div = real_t(3.75e-6);
  const int nzm = gr.nzm, nzt = gr.nzt;
  const real_t F0    = parameters_radiation.F0;
  const real_t F1    = parameters_radiation.F1;
  const real_t kappa = parameters_radiation.kappa;
  const auto zt = gr.zt, zm = gr.zm, invrs_dzt = gr.invrs_dzt;
  Real2D LWP;
  Real1D z_i =
      parameters_radiation.l_rad_above_cloud
          ? workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "simple_rad::z_i", ngrdcol)
          : Real1D{};

  // ---- Begin Code ----
  LWP = liq_water_path(workspace, ngrdcol, nzm, nzt, rho, rcm, invrs_dzt);

  if (F1 > eps) {
    Kokkos::parallel_for(
        "simple_rad::both_fluxes", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          Frad_LW(i, k) = F0 * Kokkos::exp(-kappa * LWP(i, k)) +
                          F1 * Kokkos::exp(-kappa * (LWP(i, 0) - LWP(i, k)));
        });
  } else {
    // Equivalent formula when cloud-base heating vanishes, with one exponential.
    Kokkos::parallel_for(
        "simple_rad::top_flux", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          Frad_LW(i, k) = F0 * Kokkos::exp(-kappa * LWP(i, k));
        });
  }

  if (parameters_radiation.l_rad_above_cloud) {
    // Find the rtm=8.0g/kg isopleth. Native assumes an interior crossing; its
    // commented-out invalid-profile check is not replaced with JAX clipping.
    Kokkos::parallel_for(
        "simple_rad::inversion", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
          int k_iso = 0;
          while (k_iso < nzt && rtm(i, k_iso) > real_t(8.e-3))
            ++k_iso;
          z_i(i) = lin_interpolate_two_points(real_t(8.e-3), rtm(i, k_iso), rtm(i, k_iso - 1),
                                              zt(i, k_iso), zt(i, k_iso - 1));
        });

    // Evaluate the pointwise step function where it is consumed, preserving exact ties.
    Kokkos::parallel_for(
        "simple_rad::above_cloud", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          real_t Heaviside;
          if (zm(i, k) - z_i(i) < -eps) {
            Heaviside = real_t(0);
          } else if (zm(i, k) - z_i(i) > eps) {
            Heaviside = real_t(1);
          } else {
            Heaviside = real_t(0.5);
          }
          if (Heaviside > real_t(0)) {
            Frad_LW(i, k) =
                Frad_LW(i, k) +
                rho_zm(i, k) * Cp * ls_div * Heaviside *
                    (real_t(0.25) * Kokkos::pow(zm(i, k) - z_i(i), real_t(4) / real_t(3)) +
                     z_i(i) * Kokkos::pow(zm(i, k) - z_i(i), real_t(1) / real_t(3)));
          }
        });

    // Sample the inversion height after applying the above-cloud flux term.
    if (stats.l_sample) {
      stats_update("z_inversion", z_i, stats);
    }
  }

  // The radiative heating rate is defined on thermodynamic levels.
  Kokkos::parallel_for(
      "simple_rad::heating", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        radht_LW(i, k) = (one / exner(i, k)) * (-one / (Cp * rho(i, k))) *
                         (Frad_LW(i, k + 1) - Frad_LW(i, k)) * invrs_dzt(i, k);
      });
}

void simple_rad_bomex(ScratchWorkspace& workspace, // InOut
                      const grid& gr, int ngrdcol, // In
                      Real2D radht) {              // Out
  // Compute the GCSS BOMEX radiative theta-l tendency.
  // Reference: http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html,section3.4.
  const auto zt = gr.zt;
  Kokkos::parallel_for(
      "simple_rad_bomex", Kokkos::RangePolicy<>(0, ngrdcol * gr.nzt),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        if (zt(i, k) >= real_t(0) && zt(i, k) < real_t(1500)) {
          radht(i, k) = -real_t(2.315e-5);
        } else if (zt(i, k) >= real_t(1500) && zt(i, k) < real_t(2500)) {
          radht(i, k) = -real_t(2.315e-5) + real_t(2.315e-5) * (zt(i, k) - real_t(1500)) /
                                                (real_t(2500) - real_t(1500));
        } else {
          radht(i, k) = real_t(0);
        }
      });
}

void simple_rad_lba(ScratchWorkspace& workspace,                                           // InOut
                    const grid& gr, int ngrdcol, double time_current, double time_initial, // In
                    Real2D radht,                                                          // Out
                    const simple_rad_state& simple_rad // Native module In
) {
  // Compute LBA TRMM radiation by time interpolation, then height interpolation.
  auto workspace_scope = workspace.scope();
  const auto lba_krad  = simple_rad.lba_krad;
  const auto lba_zrad  = simple_rad.lba_zrad;
  Real1D radhtz        = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "simple_rad_lba::radhtz", lba_nzrad);
  real_t a;
  const real_t time = real_t(time_current - time_initial);
  int i1, i2;

  if (time <= real_t(600)) {
    Kokkos::deep_copy(radhtz, Kokkos::subview(lba_krad, Kokkos::ALL, 0));
  } else if (time >= real_t(lba_ntimes) * real_t(600)) {
    Kokkos::deep_copy(radhtz, Kokkos::subview(lba_krad, Kokkos::ALL, lba_ntimes - 1));
  } else {
    i1 = 0;
    while (i1 < lba_ntimes - 1) {
      i2 = i1 + 1;
      if (time >= real_t(600) * real_t(i1 + 1) && time < real_t(600) * real_t(i2 + 1)) {
        a               = (time - real_t(600) * real_t(i1 + 1)) /
                          (real_t(600) * real_t(i2 + 1) - real_t(600) * real_t(i1 + 1));
        const int lower = i1, upper = i2;
        Kokkos::parallel_for(
            "simple_rad_lba::time", Kokkos::RangePolicy<>(0, lba_nzrad),
            KOKKOS_LAMBDA(const int k) {
              radhtz(k) = linear_interp_factor(a, lba_krad(k, upper), lba_krad(k, lower));
            });
        i1 = lba_ntimes - 1;
      }
      // Native sets this even after the preceding assignment to lba_ntimes.
      i1 = i2;
    }
  }

  // Preserve the source's per-column calls to the existing height interpolator.
  for (int i = 0; i < ngrdcol; ++i) {
    zlinterp_fnc(gr.nzt, lba_nzrad, Kokkos::subview(gr.zt, i, Kokkos::ALL), lba_zrad, radhtz,
                 Kokkos::subview(radht, i, Kokkos::ALL));
  }
}

void simple_rad_lba_init(int iunit, const std::string& file_path, // In
                         simple_rad_state& simple_rad) {          // Native module Out
  // Read the Native33-height/36-time radiation tables through thin file I/O.
  constexpr int per_line = 5;
  simple_rad.lba_zrad =
      Real1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "simple_rad::lba_zrad"), lba_nzrad);
  simple_rad.lba_krad =
      Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "simple_rad::lba_krad"), lba_nzrad,
             lba_ntimes);
  file_read_1d(iunit, file_path + "lba_heights.dat", lba_nzrad, per_line, simple_rad.lba_zrad);
  file_read_2d(iunit, file_path + "lba_rad.dat", lba_nzrad, lba_ntimes, per_line,
               simple_rad.lba_krad);
}

namespace {
Real2D liq_water_path(ScratchWorkspace& workspace,                // InOut
                      int ngrdcol, int nzm, int nzt,              // In
                      Real2D rho, Real2D rcm, Real2D invrs_dzt) { // In
  // Integrate the liquid water path [kg/m^2] downwards from a zero top boundary.
  // Momentum levels lie between thermodynamic rho/rcm levels on the CLUBB grid.
  Real2D liq_water_path =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "liq_water_path", ngrdcol, nzm);
  Kokkos::parallel_for(
      "liq_water_path", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        liq_water_path(i, nzm - 1) = real_t(0);
        for (int k = nzm - 2; k >= 0; --k) {
          liq_water_path(i, k) = liq_water_path(i, k + 1) + rcm(i, k) * rho(i, k) / invrs_dzt(i, k);
        }
      });
  return liq_water_path;
}
} // namespace
} // namespace clubb_kokkos
