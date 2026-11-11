// Interpret IOP observations as CLUBB initial profiles and prescribed forcings.
#include "scratch_workspace.hpp"
#include "readiopdata_module.hpp"
#include "readiopdata_core/readiopdata_core_module.hpp"
#include "Input_fields/hydrostatic_module.hpp"
#include "CLUBB_core/calc_pressure.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/saturation.hpp"
#include <stdexcept>

namespace clubb_kokkos {

void readiopdata_dims(ScratchWorkspace& workspace, real_t& lat_val, real_t& lon_val, // Out
                      readiopdata_state& iop) // Native module InOut
{
  // Load dimension/time/pressure state and return the first latitude/longitude.
  auto workspace_scope = workspace.scope();
  IopReal1D lat_in, lon_in;
  int bdate, nlat, nlon;
  readiopdata_core_dims(workspace, iop.masterproc, iop.iopfilepath,
                        iop.scamiop_from_global_cam,               // In
                        iop.ntime, iop.nlev, nlat, nlon, iop.tsec, // Out
                        bdate, iop.dplevs, lat_in, lon_in,         // Out
                        iop);                                      // Native module InOut
  const auto lat_host =
      to_host_space(workspace, Kokkos::subview(lat_in, Kokkos::pair<int, int>(0, 1)));
  const auto lon_host =
      to_host_space(workspace, Kokkos::subview(lon_in, Kokkos::pair<int, int>(0, 1)));
  lat_val    = real_t(lat_host(0));
  lon_val    = real_t(lon_host(0));
  iop.l_dims = true;
  lat_in     = IopReal1D{};
  lon_in     = IopReal1D{};
}

void readiopdata_snd(ScratchWorkspace& workspace,                                         // InOut
                     real_t zm_init, int saturation_formula,                              // In
                     int& nlevels, Real1D u, Real1D v, Real1D thlm, Real1D rtm, Real1D p, // Out
                     real_t& psfc, Real1D z, Real1D ug, Real1D vg, Real1D subs,
                     std::string& subs_type, // Out
                     readiopdata_state& iop) // Native module InOut
{
  // Return the first IOP timestep on its own pressure grid, then calculate
  // altitude by the existing inverse hydrostatic integration.
  auto workspace_scope = workspace.scope();
  if (!iop.l_dims)
    throw std::runtime_error("readiopdata_dims was not called before readiopdata_snd");
  IopReal2D usnd, vsnd, tsnd, qsnd, psnd, zsnd, ugls, vgls, wgls;
  bool wgls_holds_omega, wgls_holds_w, have_w;
  readiopdata_core_snd(workspace, iop.masterproc, iop.iopfilepath,                       // In
                       iop.scamiop_from_global_cam, iop.ntime, iop.nlev,                 // In
                       iop.dplevs, float(Rd), float(Cp), float(Lv / Cp), float(Ls / Cp), // In
                       zsnd, usnd, vsnd, tsnd, qsnd, psnd, ugls, vgls, wgls,             // Out
                       &iop.get_add_surface_data, &iop.have_omega, &wgls_holds_omega,    // Out
                       &have_w, &wgls_holds_w, subs_type, // UWM_MISC Out
                       iop);                              // Native module InOut
  if (iop.get_add_surface_data)
    nlevels = iop.nlev + 1;
  else
    nlevels = iop.nlev;
  Real1D theta =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "readiopdata_snd::theta", nlevels);
  Real1D thlm_temporary = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                       "readiopdata_snd::thlm_temporary", nlevels);
  // Use only the first timestep; conversion from Native real4 occurs here.
  Kokkos::parallel_for(
      "readiopdata_snd::first_time", nlevels, KOKKOS_LAMBDA(int k) {
        u(k)     = usnd(k, 0);
        v(k)     = vsnd(k, 0);
        theta(k) = tsnd(k, 0);
        rtm(k)   = qsnd(k, 0);
        z(k)     = zsnd(k, 0);
        ug(k)    = ugls(k, 0);
        vg(k)    = vgls(k, 0);
        subs(k)  = wgls(k, 0);
      });
  const auto surface_pressure =
      to_host_space(workspace, Kokkos::subview(psnd, 0, Kokkos::pair<int, int>(0, 1)));
  psfc = real_t(surface_pressure(0)) * pascal_per_mb;
  Kokkos::parallel_for(
      "readiopdata_snd::pressure", nlevels,
      KOKKOS_LAMBDA(int k) { p(k) = psnd(k, 0) * pascal_per_mb; });
  Real1D exner =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "readiopdata_snd::exner", nlevels);
  // Native (nlevels,1) storage is sequence-associated to (1,nlevels) in
  // calculate_thvm. Use that existing column-first API shape explicitly.
  Real2D rcm =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "readiopdata_snd::rcm", 1, nlevels);
  Real2D thvm = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "readiopdata_snd::thvm",
                                             1, nlevels);
  Kokkos::parallel_for(
      "readiopdata_snd::exner", nlevels,
      KOKKOS_LAMBDA(int k) { exner(k) = Kokkos::pow(p(k) / p0, kappa); });
  Kokkos::parallel_for(
      "readiopdata_snd::mixing_ratio", nlevels,
      KOKKOS_LAMBDA(int k) { rtm(k) = rtm(k) / g_per_kg; });
  Kokkos::parallel_for(
      "readiopdata_snd::theta", nlevels, KOKKOS_LAMBDA(int k) { thlm(k) = theta(k); });
  // Initial cloud water from potential temperature, Exner and total water.
  const real_t zero_threshold = clubb_kokkos::zero_threshold;
  Kokkos::parallel_for(
      "readiopdata_snd::cloud_water", nlevels, KOKKOS_LAMBDA(int k) {
        rcm(0, k) =
            Kokkos::max(rtm(k) - sat_mixrat_liq_api(p(k), theta(k) * exner(k), saturation_formula),
                        zero_threshold);
      });
  // theta_l = theta - [Lv/(Cp*Exner)]*rcm, preserving source grouping.
  Kokkos::parallel_for(
      "readiopdata_snd::theta_l", nlevels,
      KOKKOS_LAMBDA(int k) { thlm_temporary(k) = theta(k) - Lv / (Cp * exner(k)) * rcm(0, k); });
  // Materialize each Native RESHAPE/array expression once at the actual API
  // boundary; valid prefixes retain sequence association without static tails.
  Real2D thlm_2d       = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_snd::thlm_2d", 1, nlevels);
  Real2D rtm_2d        = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_snd::rtm_2d", 1, nlevels);
  Real2D exner_2d      = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_snd::exner_2d", 1, nlevels);
  Real2D thv_reference = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_snd::thv_reference", 1, nlevels);
  Kokkos::parallel_for(
      "readiopdata_snd::thvm_inputs", nlevels, KOKKOS_LAMBDA(int k) {
        thlm_2d(0, k)       = thlm_temporary(k);
        rtm_2d(0, k)        = rtm(k);
        exner_2d(0, k)      = exner(k);
        thv_reference(0, k) = theta(k) * Kokkos::pow(one + ep2 * (rtm(k) - rcm(0, k)), kappa);
      });
  calculate_thvm(nlevels, 1,                                    // In
                 thlm_2d, rtm_2d, rcm, exner_2d, thv_reference, // In
                 thvm);                                         // Out
  // Find altitudes of the pressure sounding levels.
  inverse_hydrostatic(workspace, psfc, zm_init, nlevels,                       // In
                      Kokkos::subview(thvm, 0, Kokkos::ALL), exner,            // In
                      Kokkos::subview(z, Kokkos::pair<int, int>(0, nlevels))); // Out
  iop.l_snd = true;
  usnd      = IopReal2D{};
  vsnd      = IopReal2D{};
  tsnd      = IopReal2D{};
  qsnd      = IopReal2D{};
  psnd      = IopReal2D{};
  zsnd      = IopReal2D{};
  ugls      = IopReal2D{};
  vgls      = IopReal2D{};
  wgls      = IopReal2D{};
  theta     = Real1D{};
  exner     = Real1D{};
  rcm       = Real2D{};
  thvm      = Real2D{};
}

void readiopdata_sfc(ScratchWorkspace& workspace, Real1D* tsfc, Real1D* ubar, Real1D* shsfc,
                     Real1D* lhsfc,
                     Real1D* timesfc,        // Optional Out
                     readiopdata_state& iop) // Native module InOut
{
  // Read surface forcing, allocating only the original PRESENT outputs.
  auto workspace_scope = workspace.scope();
  if (!(iop.l_dims && iop.l_snd))
    throw std::runtime_error(
        "readiopdata_dims or readiopdata_snd was not called before readiopdata_sfc");
  IopReal1D sstsfc, tausfc, shsfc_, lhsfc_;
  readiopdata_core_sfc(workspace, iop.masterproc, iop.iopfilepath, // In
                       iop.scamiop_from_global_cam, false,         // In
                       iop.get_add_surface_data, iop.ntime,        // In
                       sstsfc, shsfc_, lhsfc_, tausfc,             // Out
                       iop);                                       // Native module InOut
  const int ntime = iop.ntime;
  if (tsfc) {
    *tsfc             = Real1D("readiopdata_sfc::tsfc", ntime);
    const auto output = *tsfc;
    Kokkos::parallel_for(
        "readiopdata_sfc::tsfc", ntime, KOKKOS_LAMBDA(int i) { output(i) = sstsfc(i); });
  }
  if (ubar) {
    *ubar             = Real1D("readiopdata_sfc::ubar", ntime);
    const auto output = *ubar;
    Kokkos::parallel_for(
        "readiopdata_sfc::ubar", ntime, KOKKOS_LAMBDA(int i) { output(i) = tausfc(i); });
  }
  if (shsfc) {
    *shsfc            = Real1D("readiopdata_sfc::shsfc", ntime);
    const auto output = *shsfc;
    Kokkos::parallel_for(
        "readiopdata_sfc::shsfc", ntime, KOKKOS_LAMBDA(int i) { output(i) = shsfc_(i); });
  }
  if (lhsfc) {
    *lhsfc            = Real1D("readiopdata_sfc::lhsfc", ntime);
    const auto output = *lhsfc;
    Kokkos::parallel_for(
        "readiopdata_sfc::lhsfc", ntime, KOKKOS_LAMBDA(int i) { output(i) = lhsfc_(i); });
  }
  if (timesfc) {
    *timesfc          = Real1D("readiopdata_sfc::timesfc", ntime);
    const auto output = *timesfc;
    const auto tsec   = iop.tsec;
    Kokkos::parallel_for(
        "readiopdata_sfc::timesfc", ntime, KOKKOS_LAMBDA(int i) { output(i) = real_t(tsec(i)); });
  }
  iop.l_sfc = true;
  shsfc_    = IopReal1D{};
  lhsfc_    = IopReal1D{};
  sstsfc    = IopReal1D{};
  tausfc    = IopReal1D{};
}

void readiopdata_frc(ScratchWorkspace& workspace,             // InOut
                     Real1D rtm, const grid& gr,              // In
                     Real1D thlm_forcing, Real1D rtm_forcing, // Out
                     readiopdata_state& iop)                  // Native module InOut
{
  // Read forcing only once and copy its first timestep, as Native does.
  auto workspace_scope = workspace.scope();
  (void)rtm;
  (void)gr; // Source's qtm_forcing automatic local is unused too.
  if (!(iop.l_dims && iop.l_snd))
    throw std::runtime_error(
        "readiopdata_dims or readiopdata_snd was not called before readiopdata_frc");
  int nlevels; // Source assigns this only on the first call; later use is undefined.
  if (!iop.l_frc) {
    readiopdata_core_frc(workspace, iop.masterproc, iop.iopfilepath, iop.get_add_surface_data, // In
                         iop.have_omega, iop.dplevs, iop.nlev, iop.ntime,                      // In
                         iop.dtls, iop.dqls, // Out
                         iop);               // Native module InOut
    if (iop.get_add_surface_data)
      nlevels = iop.nlev + 1;
    else
      nlevels = iop.nlev;
  }
  // Native also selects nlev+1 rows when its core allocated nlev rows for the
  // surface branch. Preserve the source boundary; do not pad/repair silently.
  const auto dtls = iop.dtls, dqls = iop.dqls;
  Kokkos::parallel_for(
      "readiopdata_frc::thlm", nlevels, KOKKOS_LAMBDA(int k) { thlm_forcing(k) = dtls(k, 0); });
  Kokkos::parallel_for(
      "readiopdata_frc::rtm", nlevels, KOKKOS_LAMBDA(int k) { rtm_forcing(k) = dqls(k, 0); });
  iop.l_frc = true;
}

[[noreturn]] void task_abort() {
  // Mimic the external SAM task_abort service: terminate with Native error999.
  throw std::runtime_error("999");
}
} // namespace clubb_kokkos
