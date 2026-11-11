// MPACE-A observed forcing and surface exchange; see mpace_a.hpp.
#include "scratch_workspace.hpp"
#include "CLUBB_core/file_functions.hpp"
#include "mpace_a.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "time_dependent_input.hpp"
#include <iostream>

namespace clubb_kokkos {

//----------------------------------------------------------------------
void mpace_a_tndcy(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                   const sclr_idx_type& sclr_idx, const grid& gr, double time, Real2D p_in_Pa,
                   Real2D wm_zt, Real2D wm_zm, Real2D thlm_forcing, Real2D rtm_forcing,
                   Real2D um_hoc_grid, Real2D vm_hoc_grid, Real3D sclrm_forcing,
                   Real3D edsclrm_forcing, const mpace_a_state& input) {
  // Description:
  // Interpolate the MPACE-A temperature/water and observed-wind profiles in
  // time and space, then form tendencies and scalar copies with no subsidence.
  auto workspace_scope = workspace.scope();
  (void)sclr_dim;
  (void)edsclr_dim;
  constexpr real_t p_sfc     = real_t(101000.0); // Specified surface pressure [Pa].
  constexpr int file_nlevels = mpace_a_state::file_nlevels;
  const int nzt              = gr.nzt;
  Real1D dTdt_column   = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::dTdt_column", file_nlevels);
  Real1D dqdt_column   = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::dqdt_column", file_nlevels);
  Real1D vertT_column  = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::vertT_column", file_nlevels);
  Real1D vertq_column  = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::vertq_column", file_nlevels);
  Real1D um_column     = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::um_column", file_nlevels);
  Real1D vm_column     = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::vm_column", file_nlevels);
  Real2D dqdt_hoc_grid = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::dqdt_hoc_grid", ngrdcol, nzt);
  Real2D vertT_hoc_grid = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "mpace_a_tndcy::vertT_hoc_grid", ngrdcol, nzt);
  Real2D vertq_hoc_grid = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "mpace_a_tndcy::vertq_hoc_grid", ngrdcol, nzt);
  Real2D dTdt_hoc_grid = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "mpace_a_tndcy::dTdt_hoc_grid", ngrdcol, nzt);
  int before_time = -1, after_time = -1;
  real_t ratio;

  // Select the neighboring observed times and their interpolation ratio.
  time_select(time, mpace_a_state::file_ntimes, input.file_times, before_time,
              after_time, ratio);
  if (clubb_at_least_debug_level_api(1) && (before_time == -1 || after_time == -1)) {
    std::cerr << "file_times not sorted in mpace_a_tndcy.\n";
  }
  const auto dTdt_forcing  = input.dTdt_forcing;
  const auto dqdt_forcing  = input.dqdt_forcing;
  const auto vertT_forcing = input.vertT_forcing;
  const auto vertq_forcing = input.vertq_forcing;
  const auto um_obs        = input.um_obs;
  const auto vm_obs        = input.vm_obs;
  using Policy             = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "mpace_a_tndcy::time", Policy(0, file_nlevels), KOKKOS_LAMBDA(int k) {
        dTdt_column(k) =
            linear_interp_factor(ratio, dTdt_forcing(k, after_time), dTdt_forcing(k, before_time));
        dqdt_column(k) =
            linear_interp_factor(ratio, dqdt_forcing(k, after_time), dqdt_forcing(k, before_time));
        vertT_column(k) = linear_interp_factor(ratio, vertT_forcing(k, after_time),
                                               vertT_forcing(k, before_time));
        vertq_column(k) = linear_interp_factor(ratio, vertq_forcing(k, after_time),
                                               vertq_forcing(k, before_time));
        um_column(k) = linear_interp_factor(ratio, um_obs(k, after_time), um_obs(k, before_time));
        vm_column(k) = linear_interp_factor(ratio, vm_obs(k, after_time), vm_obs(k, before_time));
      });

  // Fill each final profile directly, retaining all six source calls per column.
  for (int i = 0; i < ngrdcol; ++i) {
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 dTdt_column, Kokkos::subview(dTdt_hoc_grid, i, Kokkos::ALL));
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 dqdt_column, Kokkos::subview(dqdt_hoc_grid, i, Kokkos::ALL));
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 vertT_column, Kokkos::subview(vertT_hoc_grid, i, Kokkos::ALL));
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 vertq_column, Kokkos::subview(vertq_hoc_grid, i, Kokkos::ALL));
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 um_column, Kokkos::subview(um_hoc_grid, i, Kokkos::ALL));
    zlinterp_fnc(nzt, file_nlevels, Kokkos::subview(gr.zt, i, Kokkos::ALL), input.file_heights,
                 vm_column, Kokkos::subview(vm_hoc_grid, i, Kokkos::ALL));
  }

  // Michael Falk removed omega forcing in September 2007. All momentum-grid
  // vertical motion is zero; the old omega/interpolation block stays inactive.
  Kokkos::deep_copy(wm_zm, real_t(0.0));

  // Compute large-scale tendencies and the selected scalar copies.
  Kokkos::parallel_for(
      "mpace_a_tndcy::profiles", Policy(0, ((ngrdcol * nzt))), KOKKOS_LAMBDA(int index) {
        const int i        = ((index % ngrdcol));
        const int k        = ((index / ngrdcol));
        wm_zt(i, k)        = real_t(0.0);
        thlm_forcing(i, k) = ((dTdt_hoc_grid(i, k) + vertT_hoc_grid(i, k)) *
                              Kokkos::pow(p_sfc / p_in_Pa(i, k), Rd / Cp)) /
                             sec_per_hr;
        rtm_forcing(i, k)  = (dqdt_hoc_grid(i, k) + vertq_hoc_grid(i, k)) / g_per_kg / sec_per_hr;

        // Test scalars with theta_l and rt if desired; Native guards by index only.
        if (sclr_idx.iisclr_thl > -1)
          sclrm_forcing(i, k, sclr_idx.iisclr_thl) = thlm_forcing(i, k);
        if (sclr_idx.iisclr_rt > -1)
          sclrm_forcing(i, k, sclr_idx.iisclr_rt) = rtm_forcing(i, k);
        if (sclr_idx.iiedsclr_thl > -1)
          edsclrm_forcing(i, k, sclr_idx.iiedsclr_thl) = thlm_forcing(i, k);
        if (sclr_idx.iiedsclr_rt > -1)
          edsclrm_forcing(i, k, sclr_idx.iiedsclr_rt) = rtm_forcing(i, k);
      });
}

//----------------------------------------------------------------------
void mpace_a_sfclyr(ScratchWorkspace& workspace, int ngrdcol, double time, Real1D rho_sfc,
                    Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D ustar, const mpace_a_state& input) {
  // Description:
  // Interpolate MPACE-A observed heat fluxes and convert them to kinematic
  // heat/water fluxes using each column's surface density; prescribe ustar.
  auto workspace_scope = workspace.scope();
  int before_time = -1, after_time = -1;
  real_t ratio;
  time_select(time, mpace_a_state::file_ntimes, input.file_times, before_time,
              after_time, ratio);
  // Retain the source sorted-time sanity diagnostic(J. Fasching, June 2008).
  if (clubb_at_least_debug_level_api(1) && (before_time == -1 || after_time == -1)) {
    std::cerr << "file_times not sorted in MPACE_A\n";
  }
  const auto latent = to_host_space(workspace, input.file_latent_ht);
  const real_t latent_heat_flx =
      linear_interp_factor(ratio, latent(after_time), latent(before_time));
  const auto sensible = to_host_space(workspace, input.file_sens_ht);
  const real_t sensible_heat_flx =
      linear_interp_factor(ratio, sensible(after_time), sensible(before_time));

  // Compute heat and moisture fluxes, then prescribe friction velocity.
  using Policy = Kokkos::RangePolicy<typename Real1D::execution_space>;
  Kokkos::parallel_for(
      "mpace_a_sfclyr", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
        wpthlp_sfc(i) = sensible_heat_flx / (rho_sfc(i) * Cp);
        wprtp_sfc(i)  = latent_heat_flx / (rho_sfc(i) * Lv);
        ustar(i)      = real_t(0.25);
      });
}

//----------------------------------------------------------------
void mpace_a_init(int iunit, const std::string& file_path, mpace_a_state& input) {
  // Description:
  // Read the observed coordinate/profile/surface files used by the tendency
  // and surface routines. Arrays replace fixed Native threadprivate fields;
  // only the raw formatted reads use the existing transient I/O bridge.
  input.file_pressure = Real1D("mpace_a::file_pressure", mpace_a_state::file_nlevels);
  input.file_heights  = Real1D("mpace_a::file_heights", mpace_a_state::file_nlevels);
  input.file_times    = decltype(input.file_times)(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mpace_a::file_times"),
      mpace_a_state::file_ntimes);
  input.dTdt_forcing =
      Real2D("mpace_a::dTdt_forcing", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.dqdt_forcing =
      Real2D("mpace_a::dqdt_forcing", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.vertT_forcing =
      Real2D("mpace_a::vertT_forcing", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.vertq_forcing =
      Real2D("mpace_a::vertq_forcing", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.um_obs = Real2D("mpace_a::um_obs", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.vm_obs = Real2D("mpace_a::vm_obs", mpace_a_state::file_nlevels, mpace_a_state::file_ntimes);
  input.file_latent_ht = Real1D("mpace_a::file_latent_ht", mpace_a_state::file_ntimes);
  input.file_sens_ht   = Real1D("mpace_a::file_sens_ht", mpace_a_state::file_ntimes);

  file_read_1d(iunit, file_path + "mpace_a_press.dat",
               mpace_a_state::file_nlevels, mpace_a_state::per_line, input.file_pressure);
  file_read_1d(iunit, file_path + "mpace_a_heights.dat",
               mpace_a_state::file_nlevels, mpace_a_state::per_line, input.file_heights);
  file_read_1d(iunit, file_path + "mpace_a_times.dat",
               mpace_a_state::file_ntimes, mpace_a_state::per_line, input.file_times);
  file_read_2d(iunit, file_path + "mpace_a_dTdt.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.dTdt_forcing);
  file_read_2d(iunit, file_path + "mpace_a_dqdt_horiz.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.dqdt_forcing);
  file_read_2d(iunit, file_path + "mpace_a_verts.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.vertT_forcing);
  file_read_2d(iunit, file_path + "mpace_a_vertq.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.vertq_forcing);
  file_read_2d(iunit, file_path + "mpace_a_um_obs.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.um_obs);
  file_read_2d(iunit, file_path + "mpace_a_vm_obs.dat",
               mpace_a_state::file_nlevels, mpace_a_state::file_ntimes,
               mpace_a_state::per_line, input.vm_obs);
  file_read_1d(iunit, file_path + "mpace_a_lh.dat",
               mpace_a_state::file_ntimes, mpace_a_state::per_line, input.file_latent_ht);
  file_read_1d(iunit, file_path + "mpace_a_sh.dat",
               mpace_a_state::file_ntimes, mpace_a_state::per_line, input.file_sens_ht);
}

} // namespace clubb_kokkos
