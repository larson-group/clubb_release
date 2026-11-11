// Required transient public Native boundaries; implementations are pending.
//
// These declarations record the precise array contracts for integration into
// the existing f2k owner. No scalar bridge is called inside a device kernel.
// Native Nc_in_cloud_to_Ncnm is elemental; the temporary adapter adds nzt and
// ngrdcol and evaluates one conformable whole-array assignment in Fortran.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace f2k {
Kokkos::View<clubb_kokkos::real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>
Nc_in_cloud_to_Ncnm(int nzt, int ngrdcol, // Shape adapter only
                    clubb_kokkos::Real2D mu_chi_1, clubb_kokkos::Real2D mu_chi_2,
                    clubb_kokkos::Real2D sigma_chi_1, clubb_kokkos::Real2D sigma_chi_2,
                    clubb_kokkos::Real2D mixt_frac, clubb_kokkos::Real2D Nc_in_cloud,
                    clubb_kokkos::Real2D cloud_frac_1, clubb_kokkos::Real2D cloud_frac_2,
                    clubb_kokkos::Real2D const_Ncnp2_on_Ncnm2,
                    clubb_kokkos::Real2D const_corr_chi_Ncn); // All In; result Out
void precip_fraction(
    const clubb_kokkos::grid& gr, int nzt, int ngrdcol, int hydromet_dim, // In
    clubb_kokkos::Real3D hydromet, clubb_kokkos::Real2D cloud_frac,
    clubb_kokkos::Real2D cloud_frac_1, Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft> l_mix_rat_hm,
    Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft> l_frozen_hm, clubb_kokkos::Real1D hydromet_tol,
    clubb_kokkos::Real2D cloud_frac_2, clubb_kokkos::Real2D ice_supersat_frac,
    clubb_kokkos::Real2D ice_supersat_frac_1, clubb_kokkos::Real2D ice_supersat_frac_2,
    clubb_kokkos::Real2D mixt_frac, clubb_kokkos::Real2D clubb_params, // In
    clubb_kokkos::err_info_type& err_info,                             // InOut
    clubb_kokkos::Real2D precip_frac, clubb_kokkos::Real2D precip_frac_1,
    clubb_kokkos::Real2D precip_frac_2, clubb_kokkos::Real1D precip_frac_tol, // Out
    stats_type& stats);                                                       // InOut
} // namespace f2k
