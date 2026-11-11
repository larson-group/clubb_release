// Shortwave radiation; see the matching header for its grid and references.
#include "scratch_workspace.hpp"
#include "rad_lwsw_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"

namespace clubb_kokkos {
void sunray_sw(ScratchWorkspace& workspace,                                                 // InOut
               int ngrdcol, int nzt, Real2D rcm, Real2D rho, real_t xi_abs, Real2D dzt,     // In
               Real2D zm, Real2D zt,                                                        // In
               real_t radius, real_t A, real_t gc, real_t Fs0, real_t omega, bool l_center, // In
               Real2D Frad_SW                                                               // Out
) {
  // Description:
  // for CLEX altocumulus case
  // Written by Geert Lenderink for implementation of Shettle and
  // Weinman's formulation for radiative flux into the EUROCS
  // stratocumulus case.
  //
  // Adapted by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
  // others for COAMPS and CLUBB.
  //
  // Subroutine to compute shortwave radiative flux.
  //
  //
  // The code for sunray_sw was slightly reconstructed in order to more
  // closely follow Geert Lenderink's code for the Duynkerke et al.
  // EUROCS case.  The original formulation used in that paper comes from
  // Shettle and Weinman case.
  //
  // Tau is now computed inside this routine, as it's not needed in
  // nov11_rad.  Tau is also computed for each layer instead of as the
  // total optical depth from the top of the domain to the bottom of the
  // layer being computed.  This makes tau zero everywhere outside of
  // cloud.  F_diff and F_dir are also zero outside of cloud.
  //
  // Comments by Michael Falk, 26 January 2005.
  //
  // ADDITION TO COMMENT BY ADAM SMITH, 28 July 2006.
  // We are attempting to compare CLUBB results with the COAMPS-LES model.
  // To keep the two codes as similar as possible, this subroutine is a
  // near-duplicate of the rad_lwsw subroutine  used in COAMPS_LES.  I
  // have  modified variable declarations as needed to match CLUBB
  // standards, and sections that do not apply to CLUBB have been removed.
  //
  // For a full explanation of how this subroutine is implemented in CLUBB,
  // please see the nov11_altocu_tndcy or jun25_altocu_tndcy subroutines
  // within gcss.F.
  //
  // This version operates natively on CLUBB's bottom-up grid ordering
  // (k=1 at surface, k=nzt at top for thermo levels, k=nzm=nzt+1 at
  // top for momentum levels) and handles multiple grid columns.
  // In CLUBB's grid, thermo level k sits between momentum levels k and
  // k+1.  Optical depth (taupath) accumulates from the top of the
  // domain downward, i.e. from momentum level nzm down to level 1.

  // REFERENCES:
  // see subroutine rad_lwsw.
  //-----------------------------------------------------------------------

  // Local per-call arrays retain Native names and shapes; no mutable module data.
  auto workspace_scope = workspace.scope();
  Real2D tau =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "sunray_sw::tau", ngrdcol, nzt);
  Real2D taude =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "sunray_sw::taude", ngrdcol, nzt);
  Real1D tauc =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::tauc", ngrdcol);
  Real1D taucde =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::taucde", ngrdcol);
  Real1D taupath =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::taupath", ngrdcol);
  Real1D exmu0 =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::exmu0", ngrdcol);
  Real1D expk =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::expk", ngrdcol);
  Real1D exmk =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::exmk", ngrdcol);
  Real1D c1 = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::c1", ngrdcol);
  Real1D c2 = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "sunray_sw::c2", ngrdcol);
  real_t omegade, t1, t2, t3;
  real_t x1, x2, x3, rk, rk2, xi_abs2, rp, alpha, beta, rtt;
  real_t xp23p, xm23p, ap23b, ff, gcde;

  //  CONSTANTS/PARAMETERS
  //
  //  values added by Michael Falk and Adam Smith
  //
  // ff     : gc^2, denoted "g^2" in Duynkerke.            Unit: none
  // gcde   : Delta-Eddington transformation of gc.
  //          Notated g-prime in Duynkerke eqn.18.         Unit: None
  //
  //-----------------------------------------------------------------------
  ff   = gc * gc;
  gcde = gc / (real_t(1) + gc);

  //  Computation of tau and omega variables
  //
  //
  // tauc   : column total optical depth                   Unit: none
  // taupath: column total Delta-Eddington optical depth.  Unit: none
  // omega  : single-scattering albedo                     Unit: none
  // omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none
  // taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none
  // taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none
  //
  // Comments by Michael Falk, 26 January 2005
  //
  //-----------------------------------------------------------------------
  omegade = (real_t(1) - ff) * omega / (real_t(1) - omega * ff);

  // Compute per-layer optical depth and ascending column totals.
  Kokkos::parallel_for(
      "sunray_sw::tau", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        tauc(i) = real_t(0);
        for (int k = 0; k < nzt; ++k) {
          tau(i, k) = three_halves * rcm(i, k) * rho(i, k) * dzt(i, k) / radius / rho_lw;
          tauc(i)   = tauc(i) + tau(i, k);
        }
      });
  Kokkos::parallel_for(
      "sunray_sw::taude", Kokkos::RangePolicy<>(0, ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        taude(i, k) = (real_t(1) - omega * ff) * tau(i, k);
      });

  //  Computation of constants for radiative transfer equations
  //
  //  These variables come from Duynkerke eqn.20, which, with slight
  //  modifications, matches Shettle and Weinman between eqns.12b and 13.
  //  Duynkerke uses slightly different formulations than Shettle and
  //  Weinman:
  //
  // F0(S&W)    = F0(Duynkerke)*pi.
  // alpha(S&W) = alpha(Duynkerke)*F0(S&W).
  // beta(S&W)  = beta(Duynkerke)*F0(S&W).
  // c1(S&W)    = c1(Duynkerke)*F0(S&W).
  // c2(S&W)    = c2(Duynkerke)*F0(S&W).
  // x3(S&W)    = x3(Duynkerke)*F0(S&W).
  // y3(S&W)    = y3(Duynkerke)*F0(S&W).
  //
  // F0 is divided out of each term in several equations, and then
  // reintroduced when the actual radiative flux is computed.  The
  // computations here follow Lenderink's original sunray_sw code which
  // uses the Duynkerke formulations for these variables.
  //
  // x1     : term 1 in k equation                         Unit: none
  // x2     : term 2 in k equation                         Unit: none
  // rk     : k equation                                   Unit: none
  // rk2    : k equation squared                           Unit: none
  // x3     : term in denominator of alpha and beta        Unit: none
  // rp     : p equation                                   Unit: none
  // alpha  : alpha equation                               Unit: none
  // beta   : beta equation                                Unit: none
  //
  // The following variables are used by Lenderink to solve for
  // Duynkerke's parameters C1 and C2.  They are all dimensionless.
  //
  // rtt    : 2/3.
  // exmu0  : exponential term used in S&W eqn.14-- originally from
  //          eqn 1 (also Salby 9.35) in the source function for SW
  //          radiation
  // expk   : one of the coefficients of C2 on the left side of Shettle
  //          and Weinman eqn.14, originally from eqn.12a and eqn.12b
  // exmk   : reciprocal of expk, one of the coefficients of C1 on the
  //          left side of Shettle and Weinman eqn.14, originally from
  //          eqn.12a and eqn.12b
  // xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13
  // xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13
  // ap23b  : right side of Shettle and Weinman eqn.13
  // t1     : the other coefficient of C1 on left side of Shettle and
  //          Weinman eqn.14
  // t2     : the other coefficient of C2 on left side of Shettle and
  //          Weinman eqn.14
  // t3     : the coefficient of exmu0 on the right side of Shettle and
  //          Weinman eqn.14
  //
  // Comments by Michael Falk, 26 January 2005
  //
  //-----------------------------------------------------------------------
  x1      = real_t(1) - omegade * gcde;
  x2      = real_t(1) - omegade;
  rk      = Kokkos::sqrt(real_t(3) * x2 * x1);
  xi_abs2 = xi_abs * xi_abs;
  rk2     = rk * rk;
  x3      = real_t(4) * (real_t(1) - rk2 * xi_abs2);
  rp      = Kokkos::sqrt(real_t(3) * x2 / x1);
  alpha   = real_t(3) * omegade * xi_abs2 * (real_t(1) + gcde * x2) / x3;
  beta    = real_t(3) * omegade * xi_abs * (real_t(1) + real_t(3) * gcde * xi_abs2 * x2) / x3;

  rtt   = real_t(2) / real_t(3);
  xp23p = real_t(1) + rtt * rp;
  xm23p = real_t(1) - rtt * rp;
  ap23b = alpha + rtt * beta;
  t1    = real_t(1) - A - rtt * (real_t(1) + A) * rp;
  t2    = real_t(1) - A + rtt * (real_t(1) + A) * rp;
  t3    = (real_t(1) - A) * alpha - rtt * (real_t(1) + A) * beta + A * xi_abs;

  // Shettle and Weinman 13 and 14, adapted into Duynkerke, give two
  // equations and two unknowns which can be solved by linear combination
  // (C2) and then substitution (C1).
  //
  // exmu0, expk, exmk, taucde are per-column since they depend on the
  // column total optical depth (tauc).
  //
  // Comments by Michael Falk, 26 January 2005
  //
  //-----------------------------------------------------------------------
  Kokkos::parallel_for(
      "sunray_sw::coefficients", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        taucde(i) = (real_t(1) - omega * ff) * tauc(i);
        exmu0(i)  = Kokkos::exp(-taucde(i) / xi_abs);
        expk(i)   = Kokkos::exp(rk * taucde(i));
        exmk(i)   = real_t(1) / expk(i);
        c2(i)     = (xp23p * t3 * exmu0(i) - t1 * ap23b * exmk(i)) /
                    (xp23p * t2 * expk(i) - xm23p * t1 * exmk(i));
        c1(i)     = (ap23b - c2(i) * xm23p) / xp23p;
      });

  //  Computation of diffuse and direct shortwave flux
  //
  // F_diff is the first term in Duynkerke eqn.19.  The F0 and pi
  // constants which are divided out in the Duynkerke formulation are
  // reintroduced here.
  //
  // Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8,
  // where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,
  // where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in
  // Shettle and Weinman eqn.12b should actually be a lowercase p.
  //
  // F_dir is the second term in Duynkerke eqn.19.
  //
  // The negative sign for F_diff and F_dir is related to the definition
  // of which direction is a positive flux.
  //
  // Comments by Michael Falk, 26 January 2005
  //
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  //
  // Computation of shortwave fluxes on staggered grid
  //
  // For a full explanation see the "Computation of radiative fluxes on
  // staggered grid" section above.  For Frad_SW to be correctly computed
  // on w levels, the non-constant component of F_diff and F_dir, we
  // compute taupath on w levels as a centered difference between tau
  // values on mass levels.
  //
  // This uses CLUBB's bottom-up grid ordering (k=1 at surface):
  //
  // --------taupath-->--Frad_SW----------------    k = nzt+1  (top m level)
  //        /                   \
  // -taude----------------------               k = nzt    (top t level)
  //        \                   /
  // --------taupath-->--Frad_SW----------------    k = nzt    (m level)
  //        /                   \
  // -taude----------------------               k = nzt-1  (t level)
  //        \                   /
  // --------taupath-->--Frad_SW----------------    k = nzt-1  (m level)
  //                        ...
  // --------taupath-->--Frad_SW----------------    k = 1      (sfc m level)
  //
  // At momentum level k, the thermo level above is k and below is k-1.
  // taupath accumulates from the top (k=nzt+1) downward to the surface
  // (k=1).
  //
  // Vince Larson changed the F variables to w levels  03 Feb 2005
  // Michael Falk changed the loop to start at k=2 and then solved
  // separately for k=1 so the array didn't go out of bounds.
  //
  // This code makes the same assumption as above that dwm=dmw.
  //
  // Comments by Michael Falk, 16 February 2005.                          c
  //                                                                      c
  // ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
  //                  compatible with the use of a stretched
  //                  (or unevenly-spaced) grid, as well as with the use
  //                  of an evenly-spaced grid.  Therefore, dwm is not
  //                  necessarily equal to dmw.  Interpolation functions
  //                  are used to compute any weighted averages, rather
  //                  than using general numbers such as (1/2), which is
  //                  compatible only with an evenly-spaced grid.
  //                  Brian Griffin; May 10, 2008.
  //
  //-----------------------------------------------------------------------
  // taupath depends on the preceding vertical result; columns are independent.
  Kokkos::parallel_for(
      "sunray_sw::flux", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        real_t F_diff_k, F_dir_k;
        // Top momentum level: Native k=nzt+1.
        if (l_center) {
          taupath(i) = real_t(0.5) * taude(i, nzt - 1);
        } else {
          taupath(i) = real_t(0);
        }
        F_diff_k =
            (-real_t(4) / real_t(3)) * Fs0 *
            (rp * (c1(i) * Kokkos::exp(-rk * taupath(i)) - c2(i) * Kokkos::exp(rk * taupath(i))) -
             beta * Kokkos::exp(-taupath(i) / xi_abs));
        F_dir_k         = -Fs0 * xi_abs * Kokkos::exp(-taupath(i) / xi_abs);
        Frad_SW(i, nzt) = F_diff_k + F_dir_k;

        // Descend from Native momentum k=nzt through k=2.
        for (int k = nzt; k >= 2; --k) {
          if (l_center) {
            taupath(i) =
                taupath(i) + lin_interpolate_two_points(zm(i, k - 1), zt(i, k - 1), zt(i, k - 2),
                                                        taude(i, k - 1), taude(i, k - 2));
          } else {
            taupath(i) = taupath(i) + taude(i, k - 2);
          }
          F_diff_k =
              (-real_t(4) / real_t(3)) * Fs0 *
              (rp * (c1(i) * Kokkos::exp(-rk * taupath(i)) - c2(i) * Kokkos::exp(rk * taupath(i))) -
               beta * Kokkos::exp(-taupath(i) / xi_abs));
          F_dir_k           = -Fs0 * xi_abs * Kokkos::exp(-taupath(i) / xi_abs);
          Frad_SW(i, k - 1) = F_diff_k + F_dir_k;
        }

        // Surface: the ghost thermodynamic level repeats taude(:,1), so the
        // centered and non-centered branches both add this lowest layer.
        taupath(i) = taupath(i) + taude(i, 0);
        F_diff_k =
            (-real_t(4) / real_t(3)) * Fs0 *
            (rp * (c1(i) * Kokkos::exp(-rk * taupath(i)) - c2(i) * Kokkos::exp(rk * taupath(i))) -
             beta * Kokkos::exp(-taupath(i) / xi_abs));
        F_dir_k       = -Fs0 * xi_abs * Kokkos::exp(-taupath(i) / xi_abs);
        Frad_SW(i, 0) = F_diff_k + F_dir_k;
      });
}
} // namespace clubb_kokkos
