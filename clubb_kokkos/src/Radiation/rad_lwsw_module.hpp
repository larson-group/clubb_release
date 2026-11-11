// Compute shortwave flux on the Native bottom-up, column-first CLUBB grid.
// sunray_sw implements the Delta-Eddington/Shettle-Weinman formulation. Optical
// depth sums ascend through thermodynamic levels; flux integration descends
// from the top momentum level to the surface. l_center selects interpolation.
//REFERENCES (for rad_lwsw and sunray_sw):
//Bluestein, H. B., 1992: Synoptic-Dynamic Meteorology in
//Midlatitudes, Volume I.  Oxford Univerity Press, 431 pp.
//Duynkerke, P.G. et al., 2005: Observations and numerical simulations
//of the diurnal cycle of the EUROCS stratocumulus case.  In
//press, special EUROCS issue, Quart. Jour. Roy. Met. Soc.
//Salby, M. L., 1996: Fundamentals of Atmospheric Physics.  Academic
//Press, 627 pp.
//Shettle, E. P., and J. A. Weinman, 1970: The transfer of solar
//irradiance through inhomogenous turbid atmospheres evaluated
//by Eddington's approximation.  J. Atm. Sci., 27, 1048-1055.
//Stevens, B. et al., 2005: Evaluation of large-eddy simulations via
//observations of nocturnal marine stratocumulus.  Submitted to
//Mon. Wea. Rev.
//Wallace, J. M., and P. V. Hobbs, 1977: Atmospheric Science: An
//Introductory Survey.  Academic Press, 467 pp.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"

namespace clubb_kokkos {
void sunray_sw(ScratchWorkspace& workspace,                                                 // InOut
               int ngrdcol, int nzt, Real2D rcm, Real2D rho, real_t xi_abs, Real2D dzt,     // In
               Real2D zm, Real2D zt,                                                        // In
               real_t radius, real_t A, real_t gc, real_t Fs0, real_t omega, bool l_center, // In
               Real2D Frad_SW                                                               // Out
);

template <class Cloud, class Density, class Dz, class Zm, class Zt, class Flux>
void sunray_sw(ScratchWorkspace& workspace,                                                 // InOut
               int ngrdcol, int nzt, Cloud rcm, Density rho, real_t xi_abs, Dz dzt,         // In
               Zm zm, Zt zt,                                                                // In
               real_t radius, real_t A, real_t gc, real_t Fs0, real_t omega, bool l_center, // In
               Flux Frad_SW                                                                 // Out
) {
  // Preserve the public field surface for HostSpace and strided column windows.
  auto workspace_scope = workspace.scope();
  auto flux            = to_execution_space(workspace, Frad_SW, false);
  sunray_sw(workspace, ngrdcol, nzt, to_execution_space(workspace, rcm),
            to_execution_space(workspace, rho),                                   // In
            xi_abs, to_execution_space(workspace, dzt),                           // In
            to_execution_space(workspace, zm), to_execution_space(workspace, zt), // In
            radius, A, gc, Fs0, omega, l_center,                                  // In
            flux);                                                                // Out
  from_execution_space(workspace, Frad_SW, flux);
}
} // namespace clubb_kokkos
