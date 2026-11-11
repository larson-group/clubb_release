// SCAM IOP input interpretation, real4 science and Native calendar behavior.
#include "scratch_workspace.hpp"
#include "readiopdata_core_module.hpp"
#include "readiopdata_module.hpp"
#include "netcdf.hpp"
#include <array>
#include <vector>
#include <iostream>
namespace clubb_kokkos {
namespace {
constexpr float missing_value = -99999.f;
}
static bool isLeapYear(int iyear);

void readiopdata_core_dims(ScratchWorkspace& workspace, bool masterproc,
                           const std::string& iopfilepath,
                           bool scamiop_from_global_cam,                                // In
                           int& ntime, int& nlev, int& nlat, int& nlon, IopInt1D& tsec, // Out
                           int& bdate, IopReal1D& dplevs, IopReal1D& lat_in,
                           IopReal1D& lon_in,      // Out
                           readiopdata_state& iop) // Native module InOut
{
  // Allocate and read IOP dimensions, integer times/date, pressures and lat/lon.
  // UWM_MISC: prefer time's dimension but tsec's seconds-valued variable.
  auto workspace_scope = workspace.scope();
  iop.masterproc_      = masterproc;
  int ncid, status, varID;
  bool use_NF90_real; // Native unused helper input; bind by reference, never read.
  if (masterproc)
    std::cout << "Opening " << iopfilepath << '\n';
  status = netcdf::open(iopfilepath, ncid);
  if (status != netcdf::nc_noerr) {
    if (masterproc)
      std::cout << "ERROR(readiopdata_core_module.f90):Cant open iop dataset: " << iopfilepath
                << '\n';
    task_abort();
  }
  // Steffen Domke2019: supplied files contain both time and tsec. The dimension
  // comes from time, while tsec is preferred because time may be in day units.
  get_netcdf_dimlength(ncid, "time", ntime, status, false, iop);
  if (status != netcdf::nc_noerr)
    get_netcdf_dimlength(ncid, "tsec", ntime, status, true, iop);
  if (status != netcdf::nc_noerr) {
    if (masterproc)
      std::cout << "Error(readiopdata_core_module.f90):Could not read length of time dimension.\n";
    status = netcdf::close(ncid);
    task_abort();
  }
  status = netcdf::inq_varid(ncid, "tsec", varID);
  if (status != netcdf::nc_noerr)
    status = netcdf::inq_varid(ncid, "time", varID);
  if (status != netcdf::nc_noerr) {
    std::cout << "ERROR(readiopdata_core_module.f90):Could not get variable ID for tsec\n";
    status = netcdf::close(ncid);
    task_abort();
  }
  // Allocation failures propagate through Kokkos allocation, rather than a new
  // Native allocation-status variable. No successful science path is changed.
  tsec = IopInt1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_dims::tsec"),
                  ntime);
  Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> tsec_host =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
          Kokkos::WithoutInitializing, "readiopdata_core_dims::tsec_host", ntime);
  status = netcdf::get_var(ncid, varID, ntime, tsec_host.data());
  Kokkos::deep_copy(tsec, tsec_host);
  if (status != netcdf::nc_noerr) {
    std::cout << "ERROR(readiopdata_core_module.f90):Could not read variable tsec\n";
    status = netcdf::close(ncid);
    task_abort();
  }
  // Load the initial date, retaining the old bdate fallback.
  status = netcdf::inq_varid(ncid, "nbdate", varID);
  if (status != netcdf::nc_noerr) {
    status = netcdf::inq_varid(ncid, "bdate", varID);
    if (status != netcdf::nc_noerr) {
      std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable ID for bdate\n";
      status = netcdf::close(ncid);
      task_abort();
    }
  }
  status = netcdf::get_var(ncid, varID, bdate);
  if (status != netcdf::nc_noerr) {
    std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable bdate\n";
    status = netcdf::close(ncid);
    task_abort();
  }
  // Pressure levels are in Pa. Leave the extra surface slot zero, as Native.
  get_netcdf_dimlength(ncid, "lev", nlev, status, true, iop);
  dplevs = IopReal1D("readiopdata_core_dims::dplevs", nlev + 1);
  Kokkos::deep_copy(dplevs, 0.f);
  get_netcdf_var1d_real(workspace, ncid, "lev", dplevs, use_NF90_real, status, true, iop);
  if (scamiop_from_global_cam)
    Kokkos::parallel_for(
        "readiopdata_core_dims::cam_pressure", nlev + 1,
        KOKKOS_LAMBDA(int n) { dplevs(n) = dplevs(n) * 100.f; });
  get_netcdf_dimlength(ncid, "lat", nlat, status, true, iop);
  get_netcdf_dimlength(ncid, "lon", nlon, status, true, iop);
  lat_in = IopReal1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_dims::lat"),
                     nlat);
  lon_in = IopReal1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_dims::lon"),
                     nlon);
  get_netcdf_var1d_real(workspace, ncid, "lat", lat_in, use_NF90_real, status, false, iop);
  get_netcdf_var1d_real(workspace, ncid, "lon", lon_in, use_NF90_real, status, false, iop);
  status = netcdf::close(ncid);
}

void readiopdata_core_snd(ScratchWorkspace& workspace, bool masterproc,
                          const std::string& iopfilepath,
                          bool scamiop_from_global_cam, // In
                          int ntime, int nlev, IopReal1D dplevs, float rgas, float cp,
                          float fac_cond, float fac_sub,                                      // In
                          IopReal2D& zsnd, IopReal2D& usnd, IopReal2D& vsnd, IopReal2D& tsnd, // Out
                          IopReal2D& qsnd, IopReal2D& psnd, IopReal2D& ugls, IopReal2D& vgls,
                          IopReal2D& wgls, // Out
                          bool* get_add_surface_data, bool* have_omega,
                          bool* wgls_holds_omega,                                   // Optional Out
                          bool* have_w, bool* wgls_holds_w, std::string& subs_type, // UWM_MISC Out
                          readiopdata_state& iop) // Native module InOut
{
  // Read initial IOP conditions and large-scale soundings on their own grid.
  // John Truesdale1996/1998; Peter Blossey2006/2008 (SAM6.5/6.7).
  // Surface input is appended only if pressure exceeds every sounding level.
  // CLDLIQ/CLDICE are added to total water, and latent-temperature corrections
  // let the first model saturation adjustment produce the initial cloud.
  // Recommended inputs are Tli=Tabs-(L/Cp)ql-(Ls/Cp)qi and qtot=qv+ql+qi.
  // Native OPTIONAL flags are dereferenced unconditionally; the supported
  // wrapper always supplies all five. Absent source-invalid outputs gain no default.
  auto workspace_scope = workspace.scope();
  iop.masterproc_      = masterproc;
  int ncid, status;
  bool use_NF90_real;
  if (masterproc)
    std::cout << "Opening " << iopfilepath << "\n";
  status = netcdf::open(iopfilepath, ncid);
  if (status != netcdf::nc_noerr) {
    if (masterproc)
      std::cout << "ERROR(readiopdata_core_module.f90):Cant open iop dataset: " << iopfilepath
                << "\n";
    task_abort();
  }
  IopReal1D Tg_in   = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_snd::Tg_in", ntime);
  IopReal1D Ts_in   = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_snd::Ts_in", ntime);
  IopReal1D Ps_in   = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_snd::Ps_in", ntime);
  IopReal1D tmp_srf = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_snd::tmp_srf", ntime);
  Kokkos::deep_copy(Tg_in, missing_value);
  Kokkos::deep_copy(Ts_in, missing_value);
  Kokkos::deep_copy(Ps_in, missing_value);
  // Surface variables are optional except pressure. UWM CAM files may contain
  // an unusable all-zero Tsair; preserve Native ANY(zero) rejection.
  get_netcdf_var1d_real(workspace, ncid, "Tsair", Ts_in, use_NF90_real, status, false, iop);
  bool have_tsair = status == netcdf::nc_noerr;
  if (scamiop_from_global_cam) {
    int zeros = 0;
    Kokkos::parallel_reduce(
        "readiopdata_core_snd::Tsair_zero", ntime,
        KOKKOS_LAMBDA(int i, int& count) {
          if (Ts_in(i) == 0.f)
            count += 1;
        },
        zeros);
    if (zeros != 0)
      have_tsair = false;
  }
  get_netcdf_var1d_real(workspace, ncid, "Tg", Tg_in, use_NF90_real, status, false, iop);
  const bool have_tg = status == netcdf::nc_noerr;
  get_netcdf_var1d_real(workspace, ncid, "Ps", Ps_in, use_NF90_real, status, true, iop);
  // Ordered extrema preserve the source comparisons without a floating sum.
  IopReal1D extrema = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_snd::extrema", 2);
  Kokkos::parallel_for(
      "readiopdata_core_snd::surface_pressure", 1, KOKKOS_LAMBDA(int) {
        float minimum_pressure = Kokkos::finite_max_v<float>;
        float maximum_level    = -Kokkos::finite_max_v<float>;
        for (int i = 0; i < ntime; ++i)
          minimum_pressure = Kokkos::min(minimum_pressure, Ps_in(i));
        for (int n = 0; n < nlev; ++n)
          maximum_level = Kokkos::max(maximum_level, dplevs(n));
        extrema(0) = minimum_pressure;
        extrema(1) = maximum_level;
      });
  const auto extrema_host = to_host_space(workspace, extrema);
  int nlev_in;
  if (extrema_host(0) <= extrema_host(1)) {
    *get_add_surface_data = false;
    nlev_in               = nlev;
  } else {
    *get_add_surface_data = true;
    nlev_in               = nlev + 1;
  }
  // Convert pressure levels and surface pressure from Pa to millibar.
  IopReal1D levs = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                   "readiopdata_core_snd::levs", nlev + 1);
  Kokkos::parallel_for(
      "readiopdata_core_snd::levs", nlev, KOKKOS_LAMBDA(int n) { levs(n) = dplevs(n) / 100.f; });
  Kokkos::parallel_for(
      "readiopdata_core_snd::Ps", ntime, KOKKOS_LAMBDA(int i) { Ps_in(i) = Ps_in(i) / 100.f; });
  IopReal2D T_in  = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::T_in", ntime, nlev_in);
  IopReal2D q_in  = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::q_in", ntime, nlev_in);
  IopReal2D u_in  = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::u_in", ntime, nlev_in);
  IopReal2D v_in  = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::v_in", ntime, nlev_in);
  IopReal2D ug_in = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::ug_in", ntime, nlev_in);
  IopReal2D vg_in = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                    "readiopdata_core_snd::vg_in", ntime, nlev_in);
  IopReal2D omega_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_snd::omega_in", ntime, nlev_in);
  IopReal2D cldliq_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_snd::cldliq_in", ntime, nlev_in);
  IopReal2D cldice_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_snd::cldice_in", ntime, nlev_in);
  IopReal2D w_in = workspace_scope.view<IopReal2D>(Kokkos::WithoutInitializing,
                                                   "readiopdata_core_snd::w_in", ntime, nlev_in);
  Kokkos::deep_copy(T_in, missing_value);
  if (scamiop_from_global_cam)
    get_netcdf_var2d_real(workspace, ncid, "t", ntime, nlev, T_in, use_NF90_real, status, true,
                          iop);
  else
    get_netcdf_var2d_real(workspace, ncid, "T", ntime, nlev, T_in, use_NF90_real, status, true,
                          iop);
  Kokkos::deep_copy(q_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "q", ntime, nlev, q_in, use_NF90_real, status, true, iop);
  Kokkos::deep_copy(cldliq_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "CLDLIQ", ntime, nlev, cldliq_in, use_NF90_real, status,
                        false, iop);
  if (status == netcdf::nc_noerr) {
    // Add condensate to total water; subtract its latent temperature release.
    Kokkos::parallel_for(
        "readiopdata_core_snd::CLDLIQ", ntime * nlev, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ntime, n = cell / ntime;
          q_in(i, n) = q_in(i, n) + cldliq_in(i, n);
          T_in(i, n) = T_in(i, n) - fac_cond * cldliq_in(i, n);
        });
  }
  Kokkos::deep_copy(cldice_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "CLDICE", ntime, nlev, cldice_in, use_NF90_real, status,
                        false, iop);
  if (status == netcdf::nc_noerr) {
    // Add condensate to total water; subtract its latent temperature release.
    Kokkos::parallel_for(
        "readiopdata_core_snd::CLDICE", ntime * nlev, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ntime, n = cell / ntime;
          q_in(i, n) = q_in(i, n) + cldice_in(i, n);
          T_in(i, n) = T_in(i, n) - fac_sub * cldice_in(i, n);
        });
  }
  Kokkos::deep_copy(omega_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "omega", ntime, nlev, omega_in, use_NF90_real, status,
                        false, iop);
  *have_omega = status == netcdf::nc_noerr;
  if (!*have_omega) {
    Kokkos::deep_copy(w_in, missing_value);
    get_netcdf_var2d_real(workspace, ncid, "w", ntime, nlev, w_in, use_NF90_real, status, false,
                          iop);
    *have_w = status == netcdf::nc_noerr;
  }
  if (*have_omega)
    subs_type = "omega[Pa/s]";
  else if (*have_w)
    subs_type = "w[m/s]";
  if (scamiop_from_global_cam && (*have_omega || *have_w)) {
    // UWM ticket87: CAM supplies total advective tendencies. Disable its
    // vertical velocity to avoid imposing vertical advection a second time.
    Kokkos::deep_copy(omega_in, 0.f);
    Kokkos::deep_copy(w_in, 0.f);
    *have_omega = false;
    *have_w     = false;
  }
  Kokkos::deep_copy(u_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "u", ntime, nlev, u_in, use_NF90_real, status, true, iop);
  Kokkos::deep_copy(v_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "v", ntime, nlev, v_in, use_NF90_real, status, true, iop);
  Kokkos::deep_copy(ug_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "ug", ntime, nlev, ug_in, use_NF90_real, status, false,
                        iop);
  bool have_geostrophic_wind = status == netcdf::nc_noerr;
  Kokkos::deep_copy(vg_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "vg", ntime, nlev, vg_in, use_NF90_real, status, false,
                        iop);
  if (status != netcdf::nc_noerr)
    have_geostrophic_wind = false;
  if (*get_add_surface_data) {
    // Replace surface slots with available observations before extrapolating.
    if (have_tsair)
      Kokkos::parallel_for(
          "readiopdata_core_snd::Tsair_surface", ntime,
          KOKKOS_LAMBDA(int i) { T_in(i, nlev) = Ts_in(i); });
    get_netcdf_var1d_real(workspace, ncid, "Ptend", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr) {
      if (*have_w)
        Kokkos::parallel_for(
            "readiopdata_core_snd::Ptend_w", ntime,
            KOKKOS_LAMBDA(int i) { w_in(i, nlev) = tmp_srf(i); });
      if (*have_omega)
        Kokkos::parallel_for(
            "readiopdata_core_snd::Ptend_omega", ntime,
            KOKKOS_LAMBDA(int i) { omega_in(i, nlev) = tmp_srf(i); });
    }
    get_netcdf_var1d_real(workspace, ncid, "usrf", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_snd::usrf", ntime,
          KOKKOS_LAMBDA(int i) { u_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "vsrf", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_snd::vsrf", ntime,
          KOKKOS_LAMBDA(int i) { v_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "qsrf", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_snd::qsrf", ntime,
          KOKKOS_LAMBDA(int i) { q_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "ugsrf", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_snd::ugsrf", ntime,
          KOKKOS_LAMBDA(int i) { ug_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "vgsrf", tmp_srf, use_NF90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_snd::vgsrf", ntime,
          KOKKOS_LAMBDA(int i) { vg_in(i, nlev) = tmp_srf(i); });
  }
  if (!have_tg) {
    if (have_tsair) {
      std::cout << "In readiopdata(): Using Tsair for Tground\n";
      Kokkos::deep_copy(Tg_in, Ts_in);
    } else {
      std::cout << "In readiopdata(): Using lowest level T for Tground\n";
      Kokkos::parallel_for(
          "readiopdata_core_snd::ground_T", ntime,
          KOKKOS_LAMBDA(int i) { Tg_in(i) = T_in(i, nlev - 1); });
    }
  }
  status = netcdf::close(ncid);
  // Peter Blossey10Mar2008: input soundings and forcing keep their own grid.
  const int nsnd = ntime, nzsnd = nlev_in, nlsf = ntime, nzlsf = nlev_in;
  if (masterproc)
    std::cout << "sounding data: nsnd=" << nsnd << "  nzsnd=" << nzsnd << "\n";
  if (masterproc)
    std::cout << "forcing data: nlsf=" << nlsf << "  nzlsf=" << nzlsf << "\n";
  zsnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::zsnd"),
                   nzsnd, nsnd);
  usnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::usnd"),
                   nzsnd, nsnd);
  vsnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::vsnd"),
                   nzsnd, nsnd);
  tsnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::tsnd"),
                   nzsnd, nsnd);
  qsnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::qsnd"),
                   nzsnd, nsnd);
  psnd = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::psnd"),
                   nzsnd, nsnd);
  ugls = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::ugls"),
                   nzlsf, nlsf);
  vgls = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::vgls"),
                   nzlsf, nlsf);
  wgls = IopReal2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_snd::wgls"),
                   nzlsf, nlsf);
  // Absolute temperature -> potential temperature, in Native real4.
  Kokkos::parallel_for(
      "readiopdata_core_snd::theta", ntime * nlev, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ntime, n = cell / ntime;
        T_in(i, n) = T_in(i, n) * Kokkos::pow(1000.f / levs(n), rgas / cp);
      });
  if (*get_add_surface_data) {
    // levs' final slot is shared mutable scratch in Native. One sequential
    // time kernel preserves both source loops and their assignments exactly.
    const bool use_omega = *have_omega, use_w = *have_w;
    Kokkos::parallel_for(
        "readiopdata_core_snd::surface_theta", 1, KOKKOS_LAMBDA(int) {
          for (int i = 0; i < ntime; ++i) {
            if (T_in(i, nlev) != missing_value) {
              levs(nlev)    = Ps_in(i);
              T_in(i, nlev) = T_in(i, nlev) * Kokkos::pow(1000.f / levs(nlev), rgas / cp);
            }
          }
        });
    Kokkos::parallel_for(
        "readiopdata_core_snd::extrapolate", 1, KOKKOS_LAMBDA(int) {
          for (int i = 0; i < ntime; ++i) {
            levs(nlev)       = Ps_in(i);
            const float coef = (levs(nlev) - levs(nlev - 2)) / (levs(nlev - 1) - levs(nlev - 2));
            if (T_in(i, nlev) == missing_value)
              T_in(i, nlev) = (1.f - coef) * T_in(i, nlev - 2) + coef * T_in(i, nlev - 1);
            if (q_in(i, nlev) == missing_value)
              q_in(i, nlev) = (1.f - coef) * q_in(i, nlev - 2) + coef * q_in(i, nlev - 1);
            if (u_in(i, nlev) == missing_value)
              u_in(i, nlev) = (1.f - coef) * u_in(i, nlev - 2) + coef * u_in(i, nlev - 1);
            if (v_in(i, nlev) == missing_value)
              v_in(i, nlev) = (1.f - coef) * v_in(i, nlev - 2) + coef * v_in(i, nlev - 1);
            if (use_omega && omega_in(i, nlev) == missing_value)
              omega_in(i, nlev) =
                  (1.f - coef) * omega_in(i, nlev - 2) + coef * omega_in(i, nlev - 1);
            if (use_w && w_in(i, nlev) == missing_value)
              w_in(i, nlev) = (1.f - coef) * w_in(i, nlev - 2) + coef * w_in(i, nlev - 1);
            if (ug_in(i, nlev) == missing_value)
              ug_in(i, nlev) = (1.f - coef) * ug_in(i, nlev - 2) + coef * ug_in(i, nlev - 1);
            if (vg_in(i, nlev) == missing_value)
              vg_in(i, nlev) = (1.f - coef) * vg_in(i, nlev - 2) + coef * vg_in(i, nlev - 1);
          }
        });
  }
  Kokkos::deep_copy(zsnd, -999.f);
  Kokkos::parallel_for(
      "readiopdata_core_snd::profiles", nzsnd * nsnd, KOKKOS_LAMBDA(int cell) {
        const int n = cell % nzsnd, i = cell / nzsnd;
        if (n == 0)
          psnd(0, i) = Ps_in(i);
        else
          psnd(n, i) = levs(nzsnd - n - 1);
        tsnd(n, i) = T_in(i, nzsnd - n - 1);
        qsnd(n, i) = q_in(i, nzsnd - n - 1);
        usnd(n, i) = u_in(i, nzsnd - n - 1);
        vsnd(n, i) = v_in(i, nzsnd - n - 1);
      });
  // Match Marat's g/kg sounding convention.
  Kokkos::parallel_for(
      "readiopdata_core_snd::g_per_kg", nzsnd * nsnd, KOKKOS_LAMBDA(int cell) {
        const int n = cell % nzsnd, i = cell / nzsnd;
        qsnd(n, i) = 1.e3f * qsnd(n, i);
      });
  Kokkos::deep_copy(wgls, 0.f);
  if (*have_omega) {
    // Omega occupies wgls; the forcing consumer interprets its declared units.
    Kokkos::parallel_for(
        "readiopdata_core_snd::omega", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
          const int n = cell % nzlsf, i = cell / nzlsf;
          wgls(n, i) = omega_in(i, nzlsf - n - 1);
        });
    *wgls_holds_omega = true;
  } else if (*have_w) {
    Kokkos::parallel_for(
        "readiopdata_core_snd::w", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
          const int n = cell % nzlsf, i = cell / nzlsf;
          wgls(n, i) = w_in(i, nzlsf - n - 1);
        });
    *wgls_holds_w = true;
  } else {
    *wgls_holds_omega = false;
    *wgls_holds_w     = false;
    std::cout << "No large-scale vertical advection is present in the dataset\n";
    task_abort();
  }
  if (have_geostrophic_wind) {
    Kokkos::parallel_for(
        "readiopdata_core_snd::geostrophic", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
          const int n = cell % nzlsf, i = cell / nzlsf;
          ugls(n, i) = ug_in(i, nzlsf - n - 1);
          vgls(n, i) = vg_in(i, nzlsf - n - 1);
        });
  } else {
    // No complete geostrophic pair: use the wind sounding as Native does.
    Kokkos::parallel_for(
        "readiopdata_core_snd::geostrophic_fallback", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
          const int n = cell % nzlsf, i = cell / nzlsf;
          ugls(n, i) = u_in(i, nzlsf - n - 1);
          vgls(n, i) = v_in(i, nzlsf - n - 1);
        });
  }
  Tg_in     = IopReal1D{};
  Ps_in     = IopReal1D{};
  Ts_in     = IopReal1D{};
  T_in      = IopReal2D{};
  q_in      = IopReal2D{};
  u_in      = IopReal2D{};
  v_in      = IopReal2D{};
  ug_in     = IopReal2D{};
  vg_in     = IopReal2D{};
  omega_in  = IopReal2D{};
  cldliq_in = IopReal2D{};
  cldice_in = IopReal2D{};
  w_in      = IopReal2D{};
}

void readiopdata_core_sfc(ScratchWorkspace& workspace, bool masterproc,
                          const std::string& iopfilepath,
                          bool scamiop_from_global_cam,                           // In
                          bool SFC_FLX_FXD, bool get_add_surface_data, int ntime, // In
                          IopReal1D& sstsfc, IopReal1D& shsfc, IopReal1D& lhsfc,
                          IopReal1D& tausfc,      // Out
                          readiopdata_state& iop) // Native module InOut
{
  // Allocate and read surface temperature, heat fluxes and wind magnitude.
  // The two Native configuration inputs are retained though unused here.
  auto workspace_scope = workspace.scope();
  (void)scamiop_from_global_cam;
  (void)get_add_surface_data;
  iop.masterproc_ = masterproc;
  int ncid, status;
  bool use_nf90_real;
  if (masterproc)
    std::cout << "Opening " << iopfilepath << "\n";
  status = netcdf::open(iopfilepath, ncid);
  if (status != netcdf::nc_noerr) {
    if (masterproc)
      std::cout << "ERROR(readiopdata_core_module.f90):Cant open iop dataset: " << iopfilepath
                << "\n";
    task_abort();
  }
  IopReal1D Tg_in     = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                        "readiopdata_core_sfc::Tg_in", ntime);
  IopReal1D shf_in    = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                        "readiopdata_core_sfc::shf_in", ntime);
  IopReal1D lhf_in    = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                        "readiopdata_core_sfc::lhf_in", ntime);
  IopReal1D tmp_srf   = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                        "readiopdata_core_sfc::tmp_srf", ntime);
  IopReal1D tausrf_in = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                        "readiopdata_core_sfc::tausrf_in", ntime);
  Kokkos::deep_copy(shf_in, missing_value);
  Kokkos::deep_copy(lhf_in, missing_value);
  Kokkos::deep_copy(tausrf_in, missing_value);
  get_netcdf_var1d_real(workspace, ncid, "shflx", shf_in, use_nf90_real, status, false, iop);
  bool have_shflx = true;
  if (status != netcdf::nc_noerr) {
    // Old name: Native backward compatibility.
    get_netcdf_var1d_real(workspace, ncid, "sh", shf_in, use_nf90_real, status, false, iop);
    if (status != netcdf::nc_noerr)
      have_shflx = false;
  }
  get_netcdf_var1d_real(workspace, ncid, "lhflx", lhf_in, use_nf90_real, status, false, iop);
  bool have_lhflx = true;
  if (status != netcdf::nc_noerr) {
    // Old name: Native backward compatibility.
    get_netcdf_var1d_real(workspace, ncid, "lh", lhf_in, use_nf90_real, status, false, iop);
    if (status != netcdf::nc_noerr)
      have_lhflx = false;
  }
  if (SFC_FLX_FXD && (!have_lhflx || !have_shflx)) {
    if (masterproc) {
      std::cout << "ERROR(readiopdata_core_module.f90): If SFC_FLX_FXD is true, \n";
      std::cout << "          shflx and lhflx needed in SCAM iop netcdf file.\n";
    }
    task_abort();
  }
  get_netcdf_var1d_real(workspace, ncid, "Tg", Tg_in, use_nf90_real, status, false, iop);
  get_netcdf_var1d_real(workspace, ncid, "usrf", tmp_srf, use_nf90_real, status, false, iop);
  bool have_tausrf = false;
  if (status == netcdf::nc_noerr) {
    Kokkos::parallel_for(
        "readiopdata_core_sfc::usrf", ntime,
        KOKKOS_LAMBDA(int i) { tausrf_in(i) = tmp_srf(i) * tmp_srf(i); });
    have_tausrf = true;
  }
  get_netcdf_var1d_real(workspace, ncid, "vsrf", tmp_srf, use_nf90_real, status, false, iop);
  if (status == netcdf::nc_noerr) {
    if (have_tausrf)
      Kokkos::parallel_for(
          "readiopdata_core_sfc::vsrf", ntime,
          KOKKOS_LAMBDA(int i) { tausrf_in(i) = tausrf_in(i) + tmp_srf(i) * tmp_srf(i); });
  } else {
    have_tausrf = false;
  }
  status         = netcdf::close(ncid);
  const int nsfc = ntime, nlsf = ntime;
  (void)nlsf;
  if (masterproc)
    std::cout << "surface forcing data: nsfc=" << nsfc << "\n";
  sstsfc = IopReal1D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_sfc::sstsfc"), nsfc);
  shsfc  = IopReal1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_sfc::shsfc"),
                     nsfc);
  lhsfc  = IopReal1D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_sfc::lhsfc"),
                     nsfc);
  tausfc = IopReal1D(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "readiopdata_core_sfc::tausfc"), nsfc);
  // Set up surface forcings. Retain Native SQRT(usrf^2+vsrf^2), not a
  // reformulated hypotenuse; the source marks this stress interpretation FIXME.
  Kokkos::parallel_for(
      "readiopdata_core_sfc::outputs", ntime, KOKKOS_LAMBDA(int i) {
        sstsfc(i) = Tg_in(i);
        shsfc(i)  = shf_in(i);
        lhsfc(i)  = lhf_in(i);
        if (have_tausrf)
          tausfc(i) = Kokkos::sqrt(tausrf_in(i));
        else
          tausfc(i) = 0.f;
      });
  Tg_in   = IopReal1D{};
  shf_in  = IopReal1D{};
  lhf_in  = IopReal1D{};
  tmp_srf = IopReal1D{};
}

void readiopdata_core_frc(ScratchWorkspace& workspace, bool masterproc,
                          const std::string& iopfilepath, bool get_add_surface_data,
                          bool have_omega,                       // In
                          IopReal1D dplevs, int nlev, int ntime, // In
                          IopReal2D& dtls, IopReal2D& dqls,      // Out
                          readiopdata_state& iop)                // Native module InOut
{
  // Read large-scale tendencies and select horizontal or full3D forcing.
  auto workspace_scope = workspace.scope();
  iop.masterproc_      = masterproc;
  int ncid, status;
  bool use_nf90_real;
  const int nlev_in = get_add_surface_data ? nlev + 1 : nlev;
  IopReal1D levs    = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_frc::levs", nlev + 1);
  Kokkos::parallel_for(
      "readiopdata_core_frc::pressure", nlev,
      KOKKOS_LAMBDA(int n) { levs(n) = dplevs(n) / 100.f; });
  if (masterproc)
    std::cout << "Opening " << iopfilepath << "\n";
  status = netcdf::open(iopfilepath, ncid);
  if (status != netcdf::nc_noerr) {
    if (masterproc)
      std::cout << "ERROR(readiopdata_core_module.f90):Cant open iop dataset: " << iopfilepath
                << "\n";
    task_abort();
  }
  IopReal1D tmp_srf = workspace_scope.view<IopReal1D>(Kokkos::WithoutInitializing,
                                                      "readiopdata_core_frc::tmp_srf", ntime);
  IopReal2D divT_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::divT_in", ntime, nlev_in);
  IopReal2D divq_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::divq_in", ntime, nlev_in);
  IopReal2D divT3d_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::divT3d_in", ntime, nlev_in);
  IopReal2D divq3d_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::divq3d_in", ntime, nlev_in);
  IopReal2D vertdivT_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::vertdivT_in", ntime, nlev_in);
  IopReal2D vertdivq_in = workspace_scope.view<IopReal2D>(
      Kokkos::WithoutInitializing, "readiopdata_core_frc::vertdivq_in", ntime, nlev_in);
  Kokkos::deep_copy(divT_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "divT", ntime, nlev, divT_in, use_nf90_real, status, false,
                        iop);
  const bool have_divT = status == netcdf::nc_noerr;
  Kokkos::deep_copy(vertdivT_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "vertdivT", ntime, nlev, vertdivT_in, use_nf90_real,
                        status, false, iop);
  const bool have_vertdivT = status == netcdf::nc_noerr;
  Kokkos::deep_copy(divT3d_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "divT3d", ntime, nlev, divT3d_in, use_nf90_real, status,
                        false, iop);
  const bool have_divT3d = status == netcdf::nc_noerr;
  Kokkos::deep_copy(divq_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "divq", ntime, nlev, divq_in, use_nf90_real, status, false,
                        iop);
  const bool have_divq = status == netcdf::nc_noerr;
  Kokkos::deep_copy(vertdivq_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "vertdivq", ntime, nlev, vertdivq_in, use_nf90_real,
                        status, false, iop);
  const bool have_vertdivq = status == netcdf::nc_noerr;
  Kokkos::deep_copy(divq3d_in, missing_value);
  get_netcdf_var2d_real(workspace, ncid, "divq3d", ntime, nlev, divq3d_in, use_nf90_real, status,
                        false, iop);
  const bool have_divq3d = status == netcdf::nc_noerr;
  if (get_add_surface_data) {
    get_netcdf_var1d_real(workspace, ncid, "divTsrf", tmp_srf, use_nf90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::divTsrf", ntime,
          KOKKOS_LAMBDA(int i) { divT_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "vertdivTsrf", tmp_srf, use_nf90_real, status, false,
                          iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::vertdivTsrf", ntime,
          KOKKOS_LAMBDA(int i) { vertdivT_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "divT3dsrf", tmp_srf, use_nf90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::divT3dsrf", ntime,
          KOKKOS_LAMBDA(int i) { divT3d_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "divqsrf", tmp_srf, use_nf90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::divqsrf", ntime,
          KOKKOS_LAMBDA(int i) { divq_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "vertdivqsrf", tmp_srf, use_nf90_real, status, false,
                          iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::vertdivqsrf", ntime,
          KOKKOS_LAMBDA(int i) { vertdivq_in(i, nlev) = tmp_srf(i); });
    get_netcdf_var1d_real(workspace, ncid, "divq3dsrf", tmp_srf, use_nf90_real, status, false, iop);
    if (status == netcdf::nc_noerr)
      Kokkos::parallel_for(
          "readiopdata_core_frc::divq3dsrf", ntime,
          KOKKOS_LAMBDA(int i) { divq3d_in(i, nlev) = tmp_srf(i); });
  }
  status = netcdf::close(ncid);
  if (get_add_surface_data) {
    // Native never assigns levs(nlev+1) in this routine. Keep that source
    // undefined boundary explicit; do not invent a pressure or extrapolation.
    Kokkos::parallel_for(
        "readiopdata_core_frc::surface", ntime, KOKKOS_LAMBDA(int i) {
          const float coef = (levs(nlev) - levs(nlev - 2)) / (levs(nlev - 1) - levs(nlev - 2));
          if (divT_in(i, nlev) == missing_value)
            divT_in(i, nlev) = (1.f - coef) * divT_in(i, nlev - 2) + coef * divT_in(i, nlev - 1);
          if (vertdivT_in(i, nlev) == missing_value)
            vertdivT_in(i, nlev) =
                (1.f - coef) * vertdivT_in(i, nlev - 2) + coef * vertdivT_in(i, nlev - 1);
          if (divT3d_in(i, nlev) == missing_value)
            divT3d_in(i, nlev) =
                (1.f - coef) * divT3d_in(i, nlev - 2) + coef * divT3d_in(i, nlev - 1);
          if (divq_in(i, nlev) == missing_value)
            divq_in(i, nlev) = (1.f - coef) * divq_in(i, nlev - 2) + coef * divq_in(i, nlev - 1);
          if (vertdivq_in(i, nlev) == missing_value)
            vertdivq_in(i, nlev) =
                (1.f - coef) * vertdivq_in(i, nlev - 2) + coef * vertdivq_in(i, nlev - 1);
          if (divq3d_in(i, nlev) == missing_value)
            divq3d_in(i, nlev) =
                (1.f - coef) * divq3d_in(i, nlev - 2) + coef * divq3d_in(i, nlev - 1);
        });
  }
  const int nzlsf = nlev, nlsf = ntime;
  // Source outputs keep nlev rows even when surface input has nlev+1 rows.
  dtls = IopReal2D("readiopdata_core_frc::dtls", nzlsf, nlsf);
  dqls = IopReal2D("readiopdata_core_frc::dqls", nzlsf, nlsf);
  Kokkos::deep_copy(dtls, 0.f);
  Kokkos::deep_copy(dqls, 0.f);
  if (have_omega) {
    // Omega is present: use horizontal advection only.
    if (have_divT)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dtls_horizontal", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dtls(n, i) = divT_in(i, nzlsf - n - 1);
          });
    if (have_divq)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dqls_horizontal", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dqls(n, i) = divq_in(i, nzlsf - n - 1);
          });
  } else {
    // Without omega, prefer3D, then vertical+horizontal, then horizontal.
    if (have_divT3d)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dtls_0", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dtls(n, i) = divT3d_in(i, nzlsf - n - 1);
          });
    else if (have_vertdivT && have_divT)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dtls_1", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dtls(n, i) = divT_in(i, nzlsf - n - 1) + vertdivT_in(i, nzlsf - n - 1);
          });
    else if (have_divT)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dtls_2", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dtls(n, i) = divT_in(i, nzlsf - n - 1);
          });
    if (have_divq3d)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dqls_0", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dqls(n, i) = divq3d_in(i, nzlsf - n - 1);
          });
    else if (have_vertdivq && have_divq)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dqls_1", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dqls(n, i) = divq_in(i, nzlsf - n - 1) + vertdivq_in(i, nzlsf - n - 1);
          });
    else if (have_divq)
      Kokkos::parallel_for(
          "readiopdata_core_frc::dqls_2", nzlsf * nlsf, KOKKOS_LAMBDA(int cell) {
            const int n = cell % nzlsf, i = cell / nzlsf;
            dqls(n, i) = divq_in(i, nzlsf - n - 1);
          });
  }
  divT_in     = IopReal2D{};
  divT3d_in   = IopReal2D{};
  vertdivT_in = IopReal2D{};
  vertdivq_in = IopReal2D{};
  divq3d_in   = IopReal2D{};
  divq_in     = IopReal2D{};
}

void get_netcdf_dimlength(int NCID, const std::string& dimName,        // In
                          int& dimlength, int& status,                 // Out
                          bool required, const readiopdata_state& iop) // In, Native module In
{
  // Find a dimension and report its length, stopping only when required.
  int dimID;
  status = netcdf::inq_dimid(NCID, dimName, dimID);
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find dimension ID for "
                  << dimName << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note(readiopdata_core_module.f90): No dimension ID for " << dimName << '\n';
      return;
    }
  }
  std::array<char, netcdf::nc_max_name> name;
  status = netcdf::inquire_dimension(NCID, dimID, name.data(), name.size(), dimlength);
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find length of " << dimName
                  << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note - readiopdata_core_module.f90 : Could not find length of " << dimName
                  << '\n';
    }
  }
}

void get_netcdf_var1d_real(ScratchWorkspace& workspace, int NCID, const std::string& varName,
                           IopReal1D var,
                           const bool& use_NF90_REAL, // In, In, InOut, In
                           int& status, bool required,
                           const readiopdata_state& iop) // Out, In, Native module In
{
  auto workspace_scope = workspace.scope();
  // Read a real4 profile with all actual file dimensions and Native start=1.
  // The unused flag stays a reference: Native callers can pass an undefined
  // value that the original routine never reads; do not evaluate it by value.
  (void)use_NF90_REAL;
  int varID;
  status = netcdf::inq_varid(NCID, varName, varID);
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable ID for " << varName
                  << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note(readiopdata_core_module.f90): Optional variable " << varName
                  << " not found\n";
      return;
    }
  }
  // Native repeats the metadata inquiry. Keep that call/error-overwrite order.
  int xtype, ndims;
  std::array<int, netcdf::nc_max_var_dims> dimids;
  status = netcdf::inquire_variable(NCID, varID, xtype, ndims, dimids.data());
  std::vector<int> start(ndims), dimlength(ndims);
  status = netcdf::inquire_variable(NCID, varID, xtype, ndims, dimids.data());
  for (int n = 0; n < ndims; ++n) {
    int dimlen;
    std::array<char, netcdf::nc_max_name> name;
    status       = netcdf::inquire_dimension(NCID, dimids[n], name.data(), name.size(), dimlen);
    start[n]     = 1;
    dimlength[n] = dimlen;
  }
  auto host = to_host_space(workspace, var);
  status    = netcdf::get_var(NCID, varID, var.extent(0), ndims, start.data(), dimlength.data(),
                              host.data());
  Kokkos::deep_copy(var, host);
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable " << varName
                  << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note (readiopdata_core_module.f90) : Could not find " << varName << '\n';
    }
  }
}

void get_netcdf_var2d_real(ScratchWorkspace& workspace,                               // InOut
                           int NCID, const std::string& varName, int ntime, int nlev, // In
                           IopReal2D var, const bool& use_NF90_REAL,                  // InOut, In
                           int& status, bool required,
                           const readiopdata_state& iop) // Out, In, Native module In
{
  // Read a real4 matrix, then retain Native's explicit flatten/reshape permutation.
  // ntime/nlev and use_NF90_REAL are unused by the Native implementation.
  auto workspace_scope = workspace.scope();
  (void)ntime;
  (void)nlev;
  (void)use_NF90_REAL;
  int varID;
  status = netcdf::inq_varid(NCID, varName, varID);
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable ID for " << varName
                  << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note(readiopdata_core_module.f90): Optional variable " << varName
                  << " not found\n";
      return;
    }
  }
  int var_type, var_ndims;
  std::array<int, netcdf::nc_max_var_dims> var_dimIDs;
  status = netcdf::inquire_variable(NCID, varID, var_type, var_ndims, var_dimIDs.data());
  if (var_ndims > 4) {
    if (iop.masterproc_)
      std::cout << "ERROR - getncdata.f90: The input var" << varName << "has" << var_ndims
                << "dimensions\n";
    status = -1;
    return;
  }
  if (var_type != netcdf::nc_float && var_type != netcdf::nc_double && var_type != netcdf::nc_int) {
    if (iop.masterproc_)
      std::cout << "ERROR - getncdata.f90: The input var" << varName << "has unknown type"
                << var_type << '\n';
    status = -1;
    return;
  }
  if (status != netcdf::nc_noerr) {
    if (iop.masterproc_)
      std::cout << "ERROR - getncdata.f90:Cant get dimension IDs for" << varName << '\n';
    return;
  }
  // The source temporary has the actual destination shape, including any surface
  // extension. Preserve unwritten cells rather than initializing them silently.
  const int rows = var.extent(0), cols = var.extent(1);
  Kokkos::View<float**, Kokkos::LayoutLeft, Kokkos::HostSpace> var_tmp_host =
      workspace_scope.view<Kokkos::View<float**, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
          Kokkos::WithoutInitializing, "get_netcdf_var2d_real::var_tmp", rows, cols);
  int xtype, ndims;
  std::array<int, netcdf::nc_max_var_dims> dimids;
  status = netcdf::inquire_variable(NCID, varID, xtype, ndims, dimids.data());
  std::vector<int> start(ndims), dimlength(ndims);
  status = netcdf::inquire_variable(NCID, varID, xtype, ndims, dimids.data());
  for (int n = 0; n < ndims; ++n) {
    int dimlen;
    std::array<char, netcdf::nc_max_name> name;
    status       = netcdf::inquire_dimension(NCID, dimids[n], name.data(), name.size(), dimlen);
    start[n]     = 1;
    dimlength[n] = dimlen;
  }
  status = netcdf::get_var(NCID, varID, rows * cols, ndims, start.data(), dimlength.data(),
                           var_tmp_host.data());
  const auto var_tmp = to_execution_space(workspace, var_tmp_host);
  // Native reorders raw (a,c,e,g;b,d,f,h) into (a,b,c,d;e,f,g,h).
  // Its status check follows this full permutation even when the read failed.
  Kokkos::parallel_for(
      "get_netcdf_var2d_real::reshape", rows * cols, KOKKOS_LAMBDA(int cell) {
        const int row_index = cell % rows, col_index = cell / rows;
        const int reshape_index   = row_index * cols + col_index;
        var(row_index, col_index) = var_tmp(reshape_index % rows, reshape_index / rows);
      });
  if (status != netcdf::nc_noerr) {
    if (required) {
      if (iop.masterproc_)
        std::cout << "ERROR(readiopdata_core_module.f90):Could not find variable " << varName
                  << '\n';
      status = netcdf::close(NCID);
      task_abort();
    } else {
      if (iop.masterproc_)
        std::cout << "Note (readiopdata_core_module.f90) : Could not find " << varName << '\n';
    }
  }
}

void calcdate(int inDate, int inSecs,                       // In
              int& outDate, int& outSecs, float& outCalday) // Out
{
  // Calculate YYMMDD, seconds into the day and fractional calendar day from a
  // base date plus elapsed seconds. Jan1 midnight is calendar day1.
  // John Truesdale; leap-year/calendar-day modifications by Peter Blossey.
  // Integer division shifts decimal date digits; MOD extracts low-order digits.
  constexpr std::array<int, 12> ndm{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  constexpr std::array<int, 12> jdcon{0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int byear = inDate / 10000;
  if (byear < 100)
    byear = 1900 + byear;
  const int bmnth = (inDate % 10000) / 100;
  const int bday  = inDate % 100;
  if (bmnth < 1 || bmnth > 12) {
    std::cout << " CALCDATE: Invalid base month input:" << bmnth << '\n';
    task_abort();
  }
  if (bday < 1 || bday > ndm[bmnth - 1]) {
    std::cout << " CALCDATE: Invalid base day of base date input:" << bday << '\n';
    task_abort();
  }
  int jday = jdcon[bmnth - 1] + bday;
  if (bmnth > 2 && isLeapYear(byear))
    jday = jday + 1;
  int jsec                = (jday - 1) * 86400 + inSecs;
  const int secs_per_year = 86400 * 365;
  int iyear;
  // Native DO retains its terminal index when no EXIT occurs.
  for (iyear = byear; iyear <= byear + inSecs / secs_per_year; ++iyear) {
    int secs_this_year = secs_per_year;
    if (isLeapYear(iyear))
      secs_this_year = secs_this_year + 86400;
    if (jsec < secs_this_year)
      break;
    jsec = jsec - secs_this_year;
  }
  const int YY = iyear % 100; // Native two-digit year convention.
  outCalday    = 1.f + float(jsec) / 86400.f;
  int MM, i;
  for (i = 0; i < 12; ++i) {
    MM                  = i + 1;
    int days_this_month = ndm[i];
    if (i == 1 && isLeapYear(iyear))
      days_this_month = days_this_month + 1;
    if (jsec < 86400 * days_this_month)
      break;
    jsec = jsec - 86400 * days_this_month;
  }
  const int DD = jsec / 86400 + 1;
  outSecs      = jsec % 86400;
  outDate      = YY * 10000 + MM * 100 + DD;
}

static bool isLeapYear(int iyear) {
  // Native internal Gregorian leap-year predicate.
  return (iyear % 4 == 0 && iyear % 100 != 0) || iyear % 400 == 0;
}

} // namespace clubb_kokkos
