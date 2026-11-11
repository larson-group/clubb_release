// Soil/vegetation surface energy budget; see the matching public header.
#include "scratch_workspace.hpp"
#include "soil_vegetation.hpp"
#include "CLUBB_core/constants_clubb.hpp"

namespace clubb_kokkos {

void advance_soil_veg(ScratchWorkspace& workspace,                       // InOut
                      int ngrdcol, real_t dt, Real1D rho_sfc,            // In
                      Real1D Frad_SW_up_sfc, Real1D Frad_SW_down_sfc,    // In
                      Real1D Frad_LW_down_sfc,                           // In
                      Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D p_sfc, // In
                      stats_type& stats,                                 // InOut
                      Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K,   // InOut
                      Real1D veg_T_in_K                                  // InOut
) {
  //
  //     Description:
  //
  //     This subroutine updates the surface and soil temp, soil heat flux,
  //     while assuming that net radiation and turbulent heat fluxes are
  //     already available from another subroutine.
  //
  //     The surface temperature (sfc_soil_T_in_K) is calculated  from
  //     the surface energy budget.
  //
  // ************************************** 2
  //               *                   *
  //               *                   *
  //             --*         +         *--  1
  //               *                   *
  // *** surface ************+************* 0
  //                        sfc_soil_T_in_K              l (evel)
  //
  //
  //     Heat conduction (in a homogeneous medium) can be described by
  //     the equation:
  //
  //              d(sfc_soil_T_in_K)/dt=ks*d(d(sfc_soil_T_in_K)/dz)/dz
  //
  //     In which ks is the soil thermal diffusity.
  //     We consider a semi half infinite medium, initially at the
  //     constant temperature deep_soil_T_in_K. if we vary the surface temperature
  //     sfc_soil_T_in_K sinussodally in time we can deduce a relation between the
  //     surface temperature, the soil heat flux (soil_heat_flux)  and deep_soil_T_in_K:
  //
  //              d(sfc_soil_T_in_K)/dt=c1*soil_heat_flux - c2*(sfc_soil_T_in_K-deep_soil_T_in_K)
  //
  //     However in reality the temperature deep_soil_T_in_K also varies in time, it
  //     may be calculated from:
  //
  //              d(deep_soil_T_in_K)/dt= c3*soil_heat_flux
  //
  //     The equations given above are analogous to those used by
  //     Deardorff  (1978).
  //
  //     Reference:
  //     Duynkerke, Peter G.  "Radiation Fog: A Comparison of Model
  //     Simulation with Detailed Observations"
  //     (February 1991) _Monthly Weather Review_ Vol 119, p. 324-341.
  //
  //
  //-----------------------------------------------------------------------

  auto workspace_scope = workspace.scope();
  real_t cs, ks, rs, c1, c2, c3, d1;
  Real1D soil_heat_flux = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                       "advance_soil_veg::soil_heat_flux", ngrdcol);

  // --------------------------------- Begin Code ---------------------------------
  // Soil parameters.
  cs = real_t(2.00e3);  // Soil heat capacity [J/(kg K)].
  rs = real_t(1.00e3);  // Soil density [kg/m^3].
  ks = real_t(2.00e-7); // Soil thermal diffusivity [m^2/s].
  d1 = Kokkos::sqrt(ks * real_t(3600) * real_t(24));
  c1 = real_t(2) * Kokkos::sqrt(pi) / (rs * cs * d1);
  c2 = real_t(2) * pi / (real_t(3600) * real_t(24));
  c3 = Kokkos::sqrt(pi * real_t(2)) / (Kokkos::exp(pi / real_t(4)) * rs * cs *
                                       Kokkos::sqrt(ks * real_t(3600) * real_t(24) * real_t(365)));

  // Sample the three temperatures before applying this timestep's updates.
  if (stats.l_sample) {
    stats_update("veg_T_in_K", veg_T_in_K, stats);
    stats_update("sfc_soil_T_in_K", sfc_soil_T_in_K, stats);
    stats_update("deep_soil_T_in_K", deep_soil_T_in_K, stats);
  }

  Kokkos::parallel_for(
      "advance_soil_veg", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        // Native integer fourth power, expressed as the square of a square.
        const real_t Frad_LW_up_sfc =
            stefan_boltzmann * ((veg_T_in_K(i) * veg_T_in_K(i)) * (veg_T_in_K(i) * veg_T_in_K(i)));
        const real_t wpthep =
            wpthlp_sfc(i) + (Lv / Cp) * Kokkos::pow(p0 / p_sfc(i), kappa) * wprtp_sfc(i);

        // Net radiation minus turbulent heat flux [W/m^2].
        const real_t veg_heat_flux = Frad_LW_down_sfc(i) - Frad_LW_up_sfc -
                                     wpthep * rho_sfc(i) * Cp +
                                     (Frad_SW_down_sfc(i) - Frad_SW_up_sfc(i));

        // Duynkerke (1991),equation19,p.328 used3W/(m^2 K); Native uses10.
        soil_heat_flux(i) =
            real_t(10) * (veg_T_in_K(i) - sfc_soil_T_in_K(i)) + real_t(0.05) * Frad_SW_down_sfc(i);

        // Preserve vegetation, surface-soil, then deep-soil update order.
        veg_T_in_K(i) = veg_T_in_K(i) + dt * real_t(5.e-5) * (veg_heat_flux - soil_heat_flux(i));
        sfc_soil_T_in_K(i) =
            sfc_soil_T_in_K(i) +
            dt * (c1 * soil_heat_flux(i) - c2 * (sfc_soil_T_in_K(i) - deep_soil_T_in_K(i)));
        deep_soil_T_in_K(i) = deep_soil_T_in_K(i) + dt * c3 * soil_heat_flux(i);
      });

  if (stats.l_sample) {
    const auto soil_heat_flux_host = to_host_space(workspace, soil_heat_flux);
    for (int i = 0; i < ngrdcol; ++i) {
      // The C++ stats API uses zero-based columns; the Native scalar updates are ordered.
      stats_update("soil_heat_flux", soil_heat_flux_host(i), stats, i);
    }
  }
}

void initialize_soil_veg(int ngrdcol,                                                        // In
                         Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K, Real1D veg_T_in_K, // Out
                         soil_vegetation_state& soil_vegetation // Native module Out
) {
  // Set the GABLS3 reference temperatures and disable the scheme for most cases.
  Kokkos::parallel_for(
      "initialize_soil_veg", Kokkos::RangePolicy<>(0, ngrdcol), KOKKOS_LAMBDA(const int i) {
        deep_soil_T_in_K(i) = real_t(288.58);
        sfc_soil_T_in_K(i)  = real_t(300);
        veg_T_in_K(i)       = real_t(300);
      });
  soil_vegetation.l_soil_veg = false;
}

} // namespace clubb_kokkos
