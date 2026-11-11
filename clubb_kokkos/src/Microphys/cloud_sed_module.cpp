// Cloud-water sedimentation on CLUBB's staggered grid.
#include "scratch_workspace.hpp"
#include "Microphys/cloud_sed_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
namespace clubb_kokkos {
//=====================================================================
// NOTE:  ADDITION OF RAIN EFFECTS AND CLOUD WATER SEDIMENTATION ON
//        RTM AND THLM.
//=====================================================================
//
// Equations:  rtm = rvm + rcm;
//             thlm = thm - ( Lv / (Cp*exner) ) * rcm
//
// When water condenses, latent heat is given off and theta (thm)
// increases by a factor of ( Lv / (Cp*exner) ) * rcm(condensed).
// The opposite effect occurs with evaporation.
//
//======================================================================
//||     Effect      |  rvm  |  rcm  |  rtm  |  rrm  |  thm  |  thlm  ||
//|====================================================================|
//|| Sedimentation   |       |       |       |       |       |        ||
//|| Effects of      | stays | incr. | incr. | stays | stays | decr.  ||
//|| Cloud Water.    | same  |       |       | same  | same  |        ||
//|| sed_rcm > 0     |       |       |       |       |       |        ||
//|====================================================================|
//|| Evaporation     |       |       |       |       |       |        ||
//|| of rain to      | incr. | stays | incr. | decr. | decr. | decr.  ||
//|| water vapor.    |       | same  |       |       |       |        ||
//|| cond_rrainm < 0 |       |       |       |       |       |        ||
//|--------------------------------------------------------------------|
//|| Autoconversion  |       |       |       |       |       |        ||
//|| of cloud water  | stays | decr. | decr. | incr. | stays | incr.  ||
//|| to rain water.  | same  |       |       |       | same  |        ||
//|| auto_rrainm > 0 |       |       |       |       |       |        ||
//|--------------------------------------------------------------------|
//|| Accretion of    |       |       |       |       |       |        ||
//|| cloud water by  | stays | decr. | decr. | incr. | stays | incr.  ||
//|| rain water.     | same  |       |       |       | same  |        ||
//|| accr_rrainm > 0 |       |       |       |       |       |        ||
//======================================================================
//
// Note: In HOC, cond_rrainm will always be either negative or zero.
//
// Overall effects of rain and cloud water sedimentation:
//
// (drtm/dt)|_t  = (drtm/dt)|_0 + sed_rcm
//                 - cond_rrainm - auto_rrainm - accr_rrainm
//
// (dthlm/dt)|_t = (dthlm/dt)|_0 - ( Lv / (Cp*exner) ) * sed_rcm
//                 + ( Lv / (Cp*exner) )
//                   * ( cond_rrainm + auto_rrainm + accr_rrainm )
void cloud_drop_sed(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, Real2D rcm,
                    Real2D Ncm, Real2D rho_zm, Real2D rho, Real2D exner, real_t sigma_g,
                    stats_type& stats, Real2D rcm_mc, Real2D thlm_mc) {
  // Account for cloud droplet sedimentation in every column. The downward
  // lognormal droplet flux uses Stokes settling; its divergence changes
  // cloud water and liquid potential temperature while leaving vapor unchanged.
  // rcm/Ncm are mean cloud water [kg/kg] and number [num/kg]; rho [kg/m^3],
  // exner [-], sigma_g [-]. InOut tendency units are [kg/kg/s] and [K/s].
  // Reference: Ackerman et al. (2009), Eq.7, doi:10.1175/2008MWR2582.1.
  auto workspace_scope = workspace.scope();
  Real2D Fcsed =
      workspace_scope.view<Real2D>("cloud_drop_sed::Fcsed", ngrdcol, gr.nzm); // Flux [kg/(m^2 s)].
  Real2D sed_rcm = workspace_scope.view<Real2D>("cloud_drop_sed::sed_rcm", ngrdcol,
                                                gr.nzt); // Tendency [kg/kg/s].
  Real2D rcm_zm = zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, rcm);
  Real2D Ncm_zm = zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Ncm);
  const int nzm = gr.nzm, nzt = gr.nzt;
  // Define cloud water sedimentation flux on momentum levels.
  Kokkos::parallel_for(
      "cloud_drop_sed::flux", ngrdcol * nzm, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        if (k > 0 && k < nzm - 1 && rcm_zm(i, k) > zero && Ncm_zm(i, k) > zero) {
          Fcsed(i, k) =
              real_t(1.19E8) *
              Kokkos::pow(three / (four * pi * rho_lw * Ncm_zm(i, k) * rho_zm(i, k)), two_thirds) *
              Kokkos::pow(rho_zm(i, k) * rcm_zm(i, k), five / three) *
              Kokkos::exp(five * (Kokkos::log(sigma_g) * Kokkos::log(sigma_g)));
        } else {
          // Boundary conditions and non-cloud cells have zero flux.
          Fcsed(i, k) = zero;
        }
      });
  // Find drc/dt due to cloud water sedimentation flux on thermodynamic levels.
  // Multiply Fcsed by Lv for W/m^2; divide its divergence by rho for kg/kg/s.
  sed_rcm = ddzm(gr.nzm, gr.nzt, ngrdcol, gr, Fcsed);
  Kokkos::parallel_for(
      "cloud_drop_sed::density", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        sed_rcm(i, k) = (one / rho(i, k)) * sed_rcm(i, k);
      });
  if (stats.l_sample) {
    stats_update("sed_rcm", sed_rcm, stats);
    stats_update("Fcsed", Fcsed, stats);
  }
  // + thlm/rtm_microphysics -- cloud water sedimentation.
  // Code addition by Brian for cloud water sedimentation.
  Kokkos::parallel_for(
      "cloud_drop_sed::tendencies", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        rcm_mc(i, k)  = rcm_mc(i, k) + sed_rcm(i, k);
        thlm_mc(i, k) = thlm_mc(i, k) - (Lv / (Cp * exner(i, k))) * sed_rcm(i, k);
      });
}
} // namespace clubb_kokkos
