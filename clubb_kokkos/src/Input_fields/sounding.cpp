// Sounding initialization, interpolation and extended profile preparation.
#include "scratch_workspace.hpp"
#include "CLUBB_core/error_code.hpp"
#include "sounding.hpp"
#include "input_reader.hpp"
#include "input_interpret.hpp"
#include "input_names.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "namelist.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace clubb_kokkos {
namespace {
constexpr int n_snd_var = 8;
// Fortran TRIM at string/file boundaries; leading blanks remain significant.
std::string sounding_trim(std::string text) {
  while (!text.empty() && text.back() == ' ')
    text.pop_back();
  return text;
}
} // namespace

static void read_sounding_file(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
                               int nmaxsnd, real_t p_sfc,              // In
                               int saturation_formula, real_t zm_init, // In
                               int& nlevels, Real1D z, Real1D theta, Real1D rt, Real1D u,
                               Real1D v,                                          // Out
                               Real1D ug, Real1D vg, Real1D p_in_Pa, Real1D subs, // Out
                               std::string& temperature_type, std::string& altitude_type,
                               std::string& subs_type,                  // Out
                               std::vector<one_dim_read_var>& retVars); // Out
static void read_sclr_sounding_file(int sclr_dim, const sclr_idx_type& sclr_idx, int iunit,
                                    const std::string& runtype,              // In
                                    int nmaxsnd,                             // In
                                    Real2D sclr_snd,                         // InOut
                                    std::vector<one_dim_read_var>& retVars); // Out
static void read_edsclr_sounding_file(int edsclr_dim, const sclr_idx_type& sclr_idx, int iunit,
                                      const std::string& runtype, // In
                                      int nmaxsnd,                // In
                                      Real2D edsclr_snd);         // InOut
static void
convert_snd2extended_atm(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
                         int n_snd_var, real_t p_sfc, real_t zm_init, // In
                         const std::vector<one_dim_read_var>& sounding_profiles,
                         int saturation_formula,                          // In
                         extended_atmosphere_state& extended_atmosphere); // Native module InOut
static void
load_extended_std_atm(int iunit,                                       // In
                      extended_atmosphere_state& extended_atmosphere); // Native module InOut

//=============================================================================
void read_sounding(ScratchWorkspace& workspace, int ngrdcol, int sclr_dim, int edsclr_dim,
                   const sclr_idx_type& sclr_idx,                                         // In
                   const grid& gr, int iunit, const std::string& runtype, real_t zm_init, // In
                   int saturation_formula,                                                // In
                   bool l_readiopdata, bool l_modify_ic_with_cubic_int,                   // In
                   real_t& p_sfc,                                                         // InOut
                   Real2D thlm, std::string& temperature_type, Real2D rtm, Real2D um, Real2D vm,
                   Real2D ugm, // Out
                   Real2D vgm, // Out
                   std::string& altitude_type, Real2D press, std::string& subs_type,
                   Real2D wm,                                                     // Out
                   Real1D rtm_sfc, Real1D thlm_sfc, Real3D sclrm, Real3D edsclrm, // Out
                   const parameters_radiation_state& parameters_radiation, // Native module In
                   extended_atmosphere_state& extended_atmosphere,         // Native module InOut
                   readiopdata_state& iop)                                 // Native module InOut
{
  // Description:
  //   Subroutine to initialize model variables from a namelist file
  // References:
  //   None
  //------------------------------------------------------------------------

  auto workspace_scope  = workspace.scope();
  constexpr int nmaxsnd = 10000;
  int nlevels;
  Real1D z = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::z",
                                          nmaxsnd); // Altitude[m]
  Real1D theta =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::theta",
                                   nmaxsnd); // Potential or liquid potential temperature[K]
  Real1D rt      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::rt",
                                                nmaxsnd); // Total water[kg/kg]
  Real1D u       = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::u",
                                                nmaxsnd); // Zonal wind[m/s]
  Real1D v       = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::v",
                                                nmaxsnd); // Meridional wind[m/s]
  Real1D ug      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::ug",
                                                nmaxsnd); // Zonal geostrophic wind[m/s]
  Real1D vg      = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::vg",
                                                nmaxsnd); // Meridional geostrophic wind[m/s]
  Real1D p_in_Pa = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                                "read_sounding::p_in_Pa", nmaxsnd); // Pressure[Pa]
  Real1D subs    = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "read_sounding::subs",
                                                nmaxsnd); // Vertical velocity [m/s or Pa/s]
  Real2D sclr_snd, edsclr_snd;
  std::vector<one_dim_read_var> sounding_retVars(n_snd_var), sclr_sounding_retVars(sclr_dim);
  temperature_type = theta_name;
  altitude_type = z_name;
  subs_type = wm_name;
  // Determine file availability before reading so missing files fail gracefully.
  const std::string case_name = sounding_trim(runtype);
  const bool l_sounding_exists =
      std::filesystem::exists("../input/case_data/" + case_name + "_sounding.in");
  const bool l_sclr_sounding_exists =
      std::filesystem::exists("../input/case_data/" + case_name + "_sclr_sounding.in");
  const bool l_edsclr_sounding_exists =
      std::filesystem::exists("../input/case_data/" + case_name + "_edsclr_sounding.in");
  // Status Message.
  if (clubb_at_least_debug_level_api(1)) {
    std::cout << "Path to sounding: " << case_name << "_sounding.in\n";
    std::cout << "File exists? " << l_sounding_exists << '\n';
    std::cout << "Path to sclr_sounding: " << case_name << "_sclr_sounding.in\n";
    std::cout << "File exists? " << l_sclr_sounding_exists << '\n';
    std::cout << "Path to edsclr_sounding: " << case_name << "_edsclr_sounding.in\n";
    std::cout << "File exists? " << l_edsclr_sounding_exists << '\n';
  }
  if (l_readiopdata) {
    readiopdata_snd(workspace, zm_init, saturation_formula, // In
                    nlevels, u, v, theta, rt, p_in_Pa,      // Out
                    p_sfc, z, ug, vg, subs, subs_type,      // Out
                    iop);                                   // Native module InOut
  } else if (l_sounding_exists) {
    // Read SAM-like <runtype>_sounding.in.
    read_sounding_file(workspace, iunit, runtype, nmaxsnd, p_sfc,  // In
                       saturation_formula, zm_init,                // In
                       nlevels, z, theta, rt, u, v,                // Out
                       ug, vg, p_in_Pa, subs,                      // Out
                       temperature_type, altitude_type, subs_type, // Out
                       sounding_retVars);                          // Out
  } else {
    std::cerr << "Cannot open " << case_name << "_sounding.in file\n";
    throw std::runtime_error("Fatal error in read_sounding");
    // Sounding namelist is no longer used. Joshua Fasching April2009.
  }
  // Read a passive scalar sounding, if enabled.
  if (sclr_dim > 0 || edsclr_dim > 0) {
    // Allocate large arrays and initialize them to zero, as in Native.
    sclr_snd = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "read_sounding::sclr_snd",
                                            nmaxsnd, sclr_max);
    edsclr_snd = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                              "read_sounding::edsclr_snd", nmaxsnd, sclr_max);
    Kokkos::deep_copy(sclr_snd, real_t(0.0));
    Kokkos::deep_copy(edsclr_snd, real_t(0.0));
    // Read SAM-like passive and eddy-passive scalar soundings.
    if (sclr_dim > 0) {
      if (l_sclr_sounding_exists)
        read_sclr_sounding_file(sclr_dim, sclr_idx, iunit, runtype, nmaxsnd, // In
                                sclr_snd,                                    // InOut
                                sclr_sounding_retVars);                      // Out
      else
        throw std::runtime_error("Cannot open <runtype>_sclr_sounding.in file");
    }
    if (edsclr_dim > 0) {
      if (l_edsclr_sounding_exists)
        read_edsclr_sounding_file(edsclr_dim, sclr_idx, iunit, runtype, nmaxsnd, // In
                                  edsclr_snd);                                   // InOut
      else
        throw std::runtime_error("Cannot open <runtype>_edsclr_sounding.in file");
    }
  }
  if (nlevels > nmaxsnd) {
    std::cerr << "Error in sounding: nlevels > nmaxsnd\n";
    std::cerr << "nlevels = " << nlevels << "\nnmaxsnd = " << nmaxsnd << '\n';
    throw std::runtime_error("STOP in read_sounding");
  }
  // Error check: if lowest above-model-surface themodynamic grid height
  // (gr%zt(1,1)) is lower than the lowest value from the input sounding,
  // then the linear interpolation scheme will fail.
  // Stage only for host diagnostics and this scalar guard; retain first-column
  // storage-level semantics rather than reinterpreting grid boundary metadata.
  const auto z_host =
      to_host_space(workspace, Kokkos::subview(z, Kokkos::pair<int, int>(0, nlevels)));
  const auto first_zt = to_host_space(workspace, Kokkos::subview(gr.zt, 0, Kokkos::ALL));
  if (first_zt(0) < z_host(0)) {
    std::cerr << "Lowest level of input sounding, z(1), must be"
              << " below the first above-model-surface thermodynamic level, gr%zt(1,1)\n";
    std::cerr << " First sounding level z(1) = " << z_host(0) << '\n';
    std::cerr << " First thermodynamic level gr%zt(1,1) = " << first_zt(0) << '\n';
    throw std::runtime_error("STOP in read_sounding");
  }
  // First sounding level should be near ground value

  // dschanen 1 May 2007
  // We have changed this for Nov. 11 and June 25, both of which
  // begin above the ground.
  // if ( abs(z(1)) > 1.e-8_core_rknd ) then
  if (false) {
    std::cerr << "First level of input sounding must be z=0\n";
    throw std::runtime_error("STOP in read_sounding");
  } else {
    const int nzt = gr.nzt;
    Kokkos::parallel_for(
        "read_sounding::first_values", ((ngrdcol * nzt)), KOKKOS_LAMBDA(int cell) {
          const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
          um(i, k) = u(0);
          vm(i, k) = v(0);
          ugm(i, k) = ug(0);
          vgm(i, k) = vg(0);
          thlm(i, k) = theta(0);
          rtm(i, k) = rt(0);
          press(i, k) = p_in_Pa(0);
          wm(i, k) = subs(0);
          if (k == 0) {
            if (sclr_dim > 0)
              for (int sclr = 0; sclr < sclr_dim; ++sclr)
                sclrm(i, 0, sclr) = sclr_snd(0, sclr);
            if (edsclr_dim > 0)
              for (int edsclr = 0; edsclr < edsclr_dim; ++edsclr)
                edsclrm(i, 0, edsclr) = edsclr_snd(0, edsclr);
          }
        });
  }
  if (clubb_at_least_debug_level_api(1)) {
    std::cout << "Reading in sounding information\n";
    // Host stream formatting replaces list-directed WRITE. Keep every source
    // field and its order visible; restore the caller's stream precision.
    const auto old_precision = std::cout.precision();
    std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
    auto print_values = [nlevels, &workspace](const char* label, Real1D profile) {
      auto print_scope = workspace.scope();
      const auto host =
          to_host_space(workspace, Kokkos::subview(profile, Kokkos::pair<int, int>(0, nlevels)));
      std::cout << label;
      for (int k = 0; k < nlevels; ++k)
        std::cout << ' ' << host(k);
      std::cout << '\n';
    };
    print_values("z = ", z);
    print_values("u = ", u);
    print_values("v = ", v);
    print_values("ug = ", ug);
    print_values("vg = ", vg);
    print_values("theta = ", theta);
    print_values("rt = ", rt);
    print_values("p_in_Pa = ", p_in_Pa);
    print_values("subs = ", subs);
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      const auto profile = to_host_space(
          workspace, Kokkos::subview(sclr_snd, Kokkos::pair<int, int>(0, nlevels), sclr));
      std::cout << "sclr_snd(" << (((sclr + 1))) << ") = ";
      for (int k = 0; k < nlevels; ++k)
        std::cout << ' ' << profile(k);
      std::cout << '\n';
    }
    for (int edsclr = 0; edsclr < edsclr_dim; ++edsclr) {
      const auto profile = to_host_space(
          workspace, Kokkos::subview(edsclr_snd, Kokkos::pair<int, int>(0, nlevels), edsclr));
      std::cout << "edsclr_snd(" << (((edsclr + 1))) << ") = ";
      for (int k = 0; k < nlevels; ++k)
        std::cout << ' ' << profile(k);
      std::cout << '\n';
    }
    std::cout.precision(old_precision);
  }
  // Use linear interpolation from two nearest prescribed grid points
  // (one above and one below) to initialize mean quantities in the model
  // Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )
  // Preserve assignments inside every Native bracket advance. Independent
  // model cells share the source profiles; only component Views are captured.
  const int nzt = gr.nzt;
  const auto gr_zt = gr.zt;
  const int iisclr_thl = sclr_idx.iisclr_thl;
  const int iisclr_rt = sclr_idx.iisclr_rt;
  const bool l_dycoms2_rf02 = case_name == "dycoms2_rf02";
  const int no_failure = std::numeric_limits<int>::max();
  auto fatal                = workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
      Kokkos::WithoutInitializing, "read_sounding::fatal", 1);
  Kokkos::deep_copy(fatal, no_failure);
  Kokkos::parallel_for(
      "read_sounding::interpolate", ((ngrdcol * nzt)), KOKKOS_LAMBDA(int cell) {
        const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
        int k_above, km1, k00, kp1, kp2;
        k_above = 0;
        while ((z(k_above) < gr_zt(i, k))) {
          k_above = (k_above + 1);
          if (k_above > nlevels - 1) {
            Kokkos::atomic_min(&fatal(0), ((i * nzt + k)));
            return;
          }
          // k_above > nlevels
          // situation w/ cubic int. (achieve better numerical solution convergence)
          if (l_modify_ic_with_cubic_int) {
            //use Steffen's monotone cubic interpolation method to obtain
            //smoothing initial condition profile for convergence test
            //note: vertical levels in sounding file need to be not too coarse
            if ((k_above == 0)) {
              km1 = k_above;
              k00 = 0;
              kp1 = 1;
              kp2 = 2;
            } else if ((k_above == 1)) {
              km1 = 0;
              kp1 = 1;
              kp2 = 2;
              k00 = 0;
            } else {
              km1 = (k_above - 2);
              kp1 = k_above;
              kp2 = (k_above + 1);
              k00 = (k_above - 1);
              //if z(k_above) reaches at the top level in sounding profile,
              //then use the nearest levels for interpolation
              if ((z(k_above) >= z(((nlevels - 1))))) {
                km1 = (nlevels - 3);
                kp1 = (nlevels - 1);
                kp2 = (nlevels - 1);
                k00 = (nlevels - 2);
              }
            }
            um(i, k)   = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), u(km1), u(k00), u(kp1), u(kp2));
            vm(i, k)   = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), v(km1), v(k00), v(kp1), v(kp2));
            ugm(i, k)  = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), ug(km1), ug(k00), ug(kp1), ug(kp2));
            vgm(i, k)  = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), vg(km1), vg(k00), vg(kp1), vg(kp2));
            thlm(i, k) = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), theta(km1), theta(k00), theta(kp1), theta(kp2));
            rtm(i, k)  = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                           z(kp2), rt(km1), rt(k00), rt(kp1), rt(kp2));
            press(i, k) =
                mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), z(kp2),
                                  p_in_Pa(km1), p_in_Pa(k00), p_in_Pa(kp1), p_in_Pa(kp2));
            wm(i, k) = mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                         z(kp2), subs(km1), subs(k00), subs(kp1), subs(kp2));
            if (!l_dycoms2_rf02) {
              //initial condition for tracers
              if ((sclr_dim > 0)) {
                for (int sclr = 0; sclr < sclr_dim; sclr += 1) {
                  sclrm(i, k, sclr) =
                      mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                        z(kp2), sclr_snd(km1, sclr), sclr_snd(k00, sclr),
                                        sclr_snd(kp1, sclr), sclr_snd(kp2, sclr));
                }
              }
              if ((edsclr_dim > 0)) {
                for (int edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
                  edsclrm(i, k, edsclr) =
                      mono_cubic_interp(gr_zt(i, k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1),
                                        z(kp2), edsclr_snd(km1, edsclr), edsclr_snd(k00, edsclr),
                                        edsclr_snd(kp1, edsclr), edsclr_snd(kp2, edsclr));
                }
              }
            } else {
              // DYCOMS2_RF02 case (use the same treatment as in regular situation w/ linear int.)
              ugm(i, k) = um(i, k);
              vgm(i, k) = vm(i, k);
              // Passive Scalars
              // Change this if they are not equal to theta_l and rt in RF02
              if ((iisclr_thl > -1)) {
                sclrm(i, k, iisclr_thl)   = thlm(i, k);
                edsclrm(i, k, iisclr_thl) = thlm(i, k);
              }
              if ((iisclr_rt > -1)) {
                sclrm(i, k, iisclr_rt)   = rtm(i, k);
                edsclrm(i, k, iisclr_rt) = rtm(i, k);
              }
            }
          } else {
            // default model setup
            // Regular situation w/ linear int.
            if (!l_dycoms2_rf02) {
              um(i, k)    = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       u(k_above), u(((k_above - 1))));
              vm(i, k)    = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       v(k_above), v(((k_above - 1))));
              ugm(i, k)   = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       ug(k_above), ug(((k_above - 1))));
              vgm(i, k)   = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       vg(k_above), vg(((k_above - 1))));
              thlm(i, k)  = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       theta(k_above), theta(((k_above - 1))));
              rtm(i, k)   = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       rt(k_above), rt(((k_above - 1))));
              press(i, k) = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       p_in_Pa(k_above), p_in_Pa(((k_above - 1))));
              wm(i, k)    = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                       subs(k_above), subs(((k_above - 1))));
              if ((sclr_dim > 0)) {
                for (int sclr = 0; sclr < sclr_dim; sclr += 1) {
                  sclrm(i, k, sclr) = lin_interpolate_two_points(
                      gr_zt(i, k), z(k_above), z(((k_above - 1))), sclr_snd(k_above, sclr),
                      sclr_snd(((k_above - 1)), sclr));
                }
              }
              if ((edsclr_dim > 0)) {
                for (int edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
                  edsclrm(i, k, edsclr) = lin_interpolate_two_points(
                      gr_zt(i, k), z(k_above), z(((k_above - 1))), edsclr_snd(k_above, edsclr),
                      edsclr_snd(((k_above - 1)), edsclr));
                }
              }
            } else {
              // DYCOMS II RF02 case
              if ((gr_zt(i, k) < real_t(795.0))) {
                // (Wyant, et al. 2007, eq 1--4)
                um(i, k) = (real_t(3.0) + ((real_t(4.3) * gr_zt(i, k)) / real_t(1000.0)));
                // Known magic number
                vm(i, k) = (-(real_t(9.0)) + ((real_t(5.6) * gr_zt(i, k)) / real_t(1000.0)));
                // Known magic number
                ugm(i, k)  = um(i, k);
                vgm(i, k)  = vm(i, k);
                thlm(i, k) = real_t(288.3);
                rtm(i, k)  = (real_t(9.45) / g_per_kg);
                // Known magic number
                // Passive Scalars
                // Change this if they are not equal to theta_l and rt in RF02
                if ((iisclr_thl > -1)) {
                  sclrm(i, k, iisclr_thl)   = thlm(i, k);
                  edsclrm(i, k, iisclr_thl) = thlm(i, k);
                }
                if ((iisclr_rt > -1)) {
                  sclrm(i, k, iisclr_rt)   = rtm(i, k);
                  edsclrm(i, k, iisclr_rt) = rtm(i, k);
                }
                press(i, k) =
                    lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                               p_in_Pa(k_above), p_in_Pa(((k_above - 1))));
                wm(i, k) = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                      subs(k_above), subs(((k_above - 1))));
              } else {
                // (Wyant, et al. 2007, eq 1--4)
                um(i, k) = (real_t(3.0) + ((real_t(4.3) * gr_zt(i, k)) / real_t(1000.0)));
                // Known magic number
                vm(i, k) = (-(real_t(9.0)) + ((real_t(5.6) * gr_zt(i, k)) / real_t(1000.0)));
                // Known magic number
                ugm(i, k)  = um(i, k);
                vgm(i, k)  = vm(i, k);
                thlm(i, k) = (real_t(295.0) + Kokkos::pow((gr_zt(i, k) - real_t(795.0)),
                                                          (real_t(1.0) / real_t(3.0))));
                // Known magic number
                rtm(i, k) =
                    ((real_t(5.0) -
                      (real_t(3.0) * (real_t(1.0) - Kokkos::exp(((real_t(795.0) - gr_zt(i, k)) /
                                                                 real_t(500.0)))))) /
                     g_per_kg);
                // Known magic number
                // Passive Scalars
                // Same as above
                if ((iisclr_thl > -1)) {
                  sclrm(i, k, iisclr_thl)   = thlm(i, k);
                  edsclrm(i, k, iisclr_thl) = thlm(i, k);
                }
                if ((iisclr_rt > -1)) {
                  sclrm(i, k, iisclr_rt)   = rtm(i, k);
                  edsclrm(i, k, iisclr_rt) = rtm(i, k);
                }
                press(i, k) =
                    lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                               p_in_Pa(k_above), p_in_Pa(((k_above - 1))));
                wm(i, k) = lin_interpolate_two_points(gr_zt(i, k), z(k_above), z(((k_above - 1))),
                                                      subs(k_above), subs(((k_above - 1))));
              }
            }
            // runtype
          }
          // l_modify_ic_with_cubic_int
        }
        // do while ( z(k) < gr%zt(i,k) )
      });
  int fatal_host;
  Kokkos::deep_copy(fatal_host, Kokkos::subview(fatal, 0));
  if (fatal_host != no_failure) {
    // The minimum is Native's i-outer/k-inner order, independent of launch order.
    const int i    = (fatal_host / nzt);
    const auto top = to_host_space(
        workspace, Kokkos::subview(gr_zt, i, Kokkos::pair<int, int>(((nzt - 1)), nzt)));
    const auto old_flags = std::cerr.flags();
    const auto old_precision = std::cerr.precision();
    std::cerr << "STOP Not enough sounding data to initialize grid:\n";
    std::cerr << std::fixed << std::setprecision(1) << "  highest sounding level" << std::setw(7)
              << z_host(((nlevels - 1))) << '\n'
              << "  should be higher than highest thermodynamic point" << std::setw(7) << top(0)
              << '\n';
    std::cerr.flags(old_flags);
    std::cerr.precision(old_precision);
    throw std::runtime_error("STOP in sounding");
  }
  // The sounding will be initialized to thermodynamic grid levels successfully
  // as long as the thermodynamic level 1 (at altitude gr%zt(1,1)) is at or above
  // the lowest sounding level (z(1)).  However, it is advantageous to know the
  // initial surface values of a few variables, as long as the sounding extends
  // to the surface, which is found at momentum level 1 (at altitude gr%zm(1,1)).
  const auto gr_zm = gr.zm;
  Kokkos::parallel_for(
      "read_sounding::surface", ngrdcol, KOKKOS_LAMBDA(int i) {
        if (gr_zm(i, 0) < z(0)) {
          // The surface (or model lower boundary) is below the lowest sounding
          // level.  Initialize the values of rtm_sfc and thlm_sfc to negative
          // values that will be overwritten later.
          rtm_sfc(i) = real_t(-999.0);
          thlm_sfc(i) = real_t(-999.0);
        } else {
          // The surface (or model lower boundary) is above the lowest sounding
          // level.  Use linear interpolation to find the values of rtm_sfc and
          // thlm_sfc.

          // Perform a binary search to find the two sounding levels that the
          // surface (gr%zm(1,1)) is found between.  The value returned (idx) is the
          // index of the closest value greater than or equal to gr%zm(1,1).
          // The C++ search and both sounding bracket slots are zero-based.
          const int idx = binary_search(nlevels, z, gr_zm(i, 0));
          // The surface is between sounding levels idx and idx-1; find rtm_sfc.
          rtm_sfc(i) =
              lin_interpolate_two_points(gr_zm(i, 0), z((((idx + 1) - 1))), z((((idx + 1) - 2))),
                                         rt((((idx + 1) - 1))), rt((((idx + 1) - 2))));
          // Use the same sounding bracket to find thlm_sfc.
          thlm_sfc(i) =
              lin_interpolate_two_points(gr_zm(i, 0), z((((idx + 1) - 1))), z((((idx + 1) - 2))),
                                         theta((((idx + 1) - 1))), theta((((idx + 1) - 2))));
        }
      });
  if (parameters_radiation.rad_scheme == "bugsrad") {
    // Prepare extended sounding for radiation.
    // The engine remains excluded; preserve these shared preparation bodies.
    if (parameters_radiation.l_use_default_std_atmosphere)
      load_extended_std_atm(iunit, extended_atmosphere);
    else
      convert_snd2extended_atm(workspace, iunit, runtype, n_snd_var, p_sfc, zm_init,
                               sounding_retVars, saturation_formula, extended_atmosphere);
  }
  // Deallocate sounding and scalar sounding profiles.  If this doesn't happen,
  // then we'll have a memory leak.
  deallocate_one_dim_vars(n_snd_var, sounding_retVars);
  deallocate_one_dim_vars(sclr_dim, sclr_sounding_retVars);
  if (sclr_snd.is_allocated())
    sclr_snd = Real2D{};
  if (edsclr_snd.is_allocated())
    edsclr_snd = Real2D{};
}

//=============================================================================
static void read_sounding_file(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
                               int nmaxsnd, real_t p_sfc,              // In
                               int saturation_formula, real_t zm_init, // In
                               int& nlevels, Real1D z, Real1D theta, Real1D rt, Real1D u,
                               Real1D v,                                          // Out
                               Real1D ug, Real1D vg, Real1D p_in_Pa, Real1D subs, // Out
                               std::string& temperature_type, std::string& altitude_type,
                               std::string& subs_type,                 // Out
                               std::vector<one_dim_read_var>& retVars) // Out
{
  //
  //  Description: This subroutine reads in a <runtype>_sounding.in file and
  //  returns the values contained in that file.

  //-----------------------------------------------------------------------

  auto workspace_scope = workspace.scope();
  read_one_dim_file(iunit, n_snd_var,                                                // In
                    "../input/case_data/" + sounding_trim(runtype) + "_sounding.in", // In
                    retVars);                                                        // Out
  fill_blanks_one_dim_vars(n_snd_var, retVars);
  read_z_profile(workspace, n_snd_var, nmaxsnd, retVars, p_sfc, zm_init, // In
                 saturation_formula,                                     // In
                 z, p_in_Pa, altitude_type);                             // Out
  read_theta_profile(workspace, n_snd_var, nmaxsnd, retVars, temperature_type, theta);
  read_x_profile(n_snd_var, nmaxsnd, rt_name, retVars, rt);
  read_x_profile(n_snd_var, nmaxsnd, um_name, retVars, u);
  read_x_profile(n_snd_var, nmaxsnd, vm_name, retVars, v);
  read_x_profile(n_snd_var, nmaxsnd, ug_name, retVars, ug);
  read_x_profile(n_snd_var, nmaxsnd, vg_name, retVars, vg);
  read_subs_profile(workspace, n_snd_var, nmaxsnd, retVars, subs_type, subs);
  nlevels = retVars[0].values.extent(0);
}

//=============================================================================
static void read_sclr_sounding_file(int sclr_dim, const sclr_idx_type& sclr_idx, int iunit,
                                    const std::string& runtype,             // In
                                    int nmaxsnd,                            // In
                                    Real2D sclr_snd,                        // InOut
                                    std::vector<one_dim_read_var>& retVars) // Out
{
  //
  // Description: This subroutine reads in a <runtype>_sclr_sounding.in file and
  //   returns the values contained in that file.

  // References:
  //   None
  //-----------------------------------------------------------------------------------------------

  (void)nmaxsnd; // Native explicit shape is represented by the caller's View.
  read_one_dim_file(iunit, sclr_dim,
                    "../input/case_data/" + sounding_trim(runtype) + "_sclr_sounding.in", retVars);
  // Native deliberately does not call fill_blanks_one_dim_vars here.
  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    const auto name = sounding_trim(retVars[sclr].name);
    if (name == CO2_name) {
      if (sclr != sclr_idx.iisclr_CO2 && sclr_idx.iisclr_CO2 > -1)
        throw std::runtime_error("iisclr_CO2 index does not match column.");
    } else if (name == rt_name) {
      if (sclr != sclr_idx.iisclr_rt && sclr_idx.iisclr_rt > -1)
        throw std::runtime_error("iisclr_rt index does not match column.");
    } else if (name == theta_name || name == thetal_name || name == temperature_name) {
      if (sclr != sclr_idx.iisclr_thl && sclr_idx.iisclr_thl > -1)
        throw std::runtime_error("iisclr_thl index does not match column.");
    }
    const int count = retVars[sclr].values.extent(0);
    auto output     = Kokkos::subview(sclr_snd, Kokkos::pair<int, int>(0, count), sclr);
    Kokkos::deep_copy(output, retVars[sclr].values);
  }
}

//=============================================================================
static void read_edsclr_sounding_file(int edsclr_dim, const sclr_idx_type& sclr_idx, int iunit,
                                      const std::string& runtype, // In
                                      int nmaxsnd,                // In
                                      Real2D edsclr_snd)          // InOut
{
  //
  //  Description: This subroutine reads in a <runtype>_edsclr_sounding.in file and
  //  returns the values contained in that file.
  //
  //-----------------------------------------------------------------------------------------------

  (void)nmaxsnd;
  std::vector<one_dim_read_var> retVars(edsclr_dim);
  read_one_dim_file(iunit, edsclr_dim,
                    "../input/case_data/" + sounding_trim(runtype) + "_edsclr_sounding.in",
                    retVars);
  // Native deliberately does not call fill_blanks_one_dim_vars here.
  for (int edsclr = 0; edsclr < edsclr_dim; ++edsclr) {
    const auto name = sounding_trim(retVars[edsclr].name);
    if (name == CO2_name) {
      if (edsclr != sclr_idx.iiedsclr_CO2 && sclr_idx.iiedsclr_CO2 > -1)
        throw std::runtime_error("iisclr_CO2 index does not match column.");
    } else if (name == rt_name) {
      if (edsclr != sclr_idx.iiedsclr_rt && sclr_idx.iiedsclr_rt > -1)
        throw std::runtime_error("iisclr_rt index does not match column.");
    } else if (name == theta_name || name == thetal_name || name == temperature_name) {
      if (edsclr != sclr_idx.iiedsclr_thl && sclr_idx.iiedsclr_thl > -1)
        throw std::runtime_error("iisclr_thl index does not match column.");
    }
    const int count = retVars[edsclr].values.extent(0);
    auto output     = Kokkos::subview(edsclr_snd, Kokkos::pair<int, int>(0, count), edsclr);
    Kokkos::deep_copy(output, retVars[edsclr].values);
  }
  deallocate_one_dim_vars(edsclr_dim, retVars);
}

//=============================================================================
void read_profile(ScratchWorkspace& workspace, const grid& gr, const std::string& fname, // In
                  Real1D x)                                                              // Out
{
  //       Description:
  //       Subroutine to initialize one generic model variable from file
  //------------------------------------------------------------------------

  auto workspace_scope  = workspace.scope();
  constexpr int nmaxsnd = 200;
  int nlevels = 0;
  // Keep Native fixed-size /profile/ storage; parsing and interpolation are C++.
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>
      z_host = workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
          "read_profile::z", nmaxsnd),
      var_host = workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>>(
          "read_profile::var", nmaxsnd);
  namelist::Group profile("profile");
  profile.bind("nlevels", nlevels);
  profile.bind("z", z_host);
  profile.bind("var", var_host);
  if (profile.read_file(fname) != 0)
    throw std::runtime_error(profile.error());
  // Native leaves an omitted nlevels undefined; reject it before indexing.
  if (nlevels < 1)
    throw std::runtime_error("read_profile requires a positive nlevels");
  if (nlevels > nmaxsnd) {
    std::cerr << "Error in sounding: nlevels > nmaxsnd\n";
    std::cerr << "nlevels = " << nlevels << "\nnmaxsnd = " << nmaxsnd << '\n';
    throw std::runtime_error("STOP in read_profile (sounding.F90)");
  }
  const auto gr_zt_host = to_host_space(workspace, Kokkos::subview(gr.zt, 0, Kokkos::ALL));
  if (gr_zt_host(0) < z_host(0)) {
    std::cerr << "Lowest level of input sounding, z(1), must be"
              << " below the first above-model-surface thermodynamic level, gr%zt(1,1)\n";
    std::cerr << " First sounding level z(1) = " << z_host(0) << '\n';
    std::cerr << " First thermodynamic level gr%zt(1,1) = " << gr_zt_host(0) << '\n';
    throw std::runtime_error("STOP in read_profile (sounding.F90)");
  }
  const Real1D z   = to_execution_space(workspace, z_host),
               var = to_execution_space(workspace, var_host);
  const int nzt = gr.nzt, no_failure = std::numeric_limits<int>::max();
  const auto gr_zt = gr.zt;
  auto fatal       = workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
      Kokkos::WithoutInitializing, "read_profile::fatal", 1);
  Kokkos::deep_copy(fatal, no_failure);
  // Use linear interpolation from two nearest prescribed grid points
  // (one above and one below) to initialize mean quantities in the model
  // Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )
  // The generic profile entry also retains assignment inside the bracket loop;
  // Native leaves x(i) unwritten if that loop is false.
  Kokkos::parallel_for(
      "read_profile::interpolate", nzt, KOKKOS_LAMBDA(int i) {
        int k = 0;
        while (z(k) < gr_zt(0, i)) {
          ++k;
          if (k > (nlevels - 1)) {
            Kokkos::atomic_min(&fatal(0), i);
            return;
          }
          x(i) =
              lin_interpolate_two_points(gr_zt(0, i), z(k), z(((k - 1))), var(k), var(((k - 1))));
        }
      });
  int fatal_host;
  Kokkos::deep_copy(fatal_host, Kokkos::subview(fatal, 0));
  if (fatal_host != no_failure) {
    const auto old_flags = std::cerr.flags();
    const auto old_precision = std::cerr.precision();
    std::cerr << "STOP Not enough sounding data to initialize grid:\n";
    std::cerr << std::fixed << std::setprecision(1) << " Highest sounding level" << std::setw(7)
              << z_host(((nlevels - 1))) << '\n'
              << "should be higher than highest thermodynamic point" << std::setw(7)
              << gr_zt_host(((nzt - 1))) << '\n';
    std::cerr.flags(old_flags);
    std::cerr.precision(old_precision);
    std::cout << " Filename: " << fname << '\n';
    throw std::runtime_error("STOP in read_profile (sounding.F90)");
  }
}

//=============================================================================
static void
convert_snd2extended_atm(ScratchWorkspace& workspace, int iunit, const std::string& runtype,
                         int n_snd_var, real_t p_sfc, real_t zm_init, // In
                         const std::vector<one_dim_read_var>& sounding_profiles,
                         int saturation_formula,                         // In
                         extended_atmosphere_state& extended_atmosphere) // Native module InOut
{
  //
  //  Description: This subroutine converts information retrieved from the
  //    sounding files of a case into a format usable for an extended atmosphere.
  //    The extended atmosphere profile is stored in module variables.
  //
  //  References:
  //    none
  //
  //-----------------------------------------------------------------------------------------------

  auto workspace_scope        = workspace.scope();
  constexpr int n_rad_scalars = 1;
  // Former module arrays remain in the shared case-owned radiation context.
  const int extended_atmos_dim = sounding_profiles[0].values.extent(0);
  extended_atmosphere.extended_atmos_dim = extended_atmos_dim;
  extended_atmosphere.extended_alt =
      Real1D("convert_snd2extended_atm::extended_alt", extended_atmos_dim);
  extended_atmosphere.extended_T_in_K =
      Real1D("convert_snd2extended_atm::extended_T_in_K", extended_atmos_dim);
  extended_atmosphere.extended_sp_hmdty =
      Real1D("convert_snd2extended_atm::extended_sp_hmdty", extended_atmos_dim);
  extended_atmosphere.extended_p_in_mb =
      Real1D("convert_snd2extended_atm::extended_p_in_mb", extended_atmos_dim);
  extended_atmosphere.extended_o3l =
      Real1D("convert_snd2extended_atm::extended_o3l", extended_atmos_dim);
  Real1D extended_p_in_Pa =
      workspace_scope.view<Real1D>("convert_snd2extended_atm::p_in_Pa", extended_atmos_dim);
  std::vector<one_dim_read_var> extended_rad_scalars(n_rad_scalars);
  std::string altitude_type, temperature_type;
  // Either convert to pressure or from pressure.
  read_z_profile(workspace, n_snd_var, extended_atmos_dim, sounding_profiles, p_sfc, zm_init, // In
                 saturation_formula,                                                          // In
                 extended_atmosphere.extended_alt, extended_p_in_Pa, altitude_type);          // Out
  if (altitude_type == z_name) {
    std::cerr << "Fatal error in convert_snd2extended_atm.\n";
    throw std::runtime_error("Feature not implemented");
  }
  // Capture component Views only, never the owning extended-atmosphere state.
  auto extended_p_in_mb = extended_atmosphere.extended_p_in_mb;
  Kokkos::parallel_for(
      "convert_snd2extended_atm::pressure", extended_atmos_dim,
      KOKKOS_LAMBDA(int i) { extended_p_in_mb(i) = extended_p_in_Pa(i) / pascal_per_mb; });
  // Convert to temperature from theta_l or theta, exactly as Native (no extra
  // cloud-water adjustment is introduced at this source boundary).
  read_theta_profile(workspace, n_snd_var, extended_atmos_dim, sounding_profiles, temperature_type,
                     extended_atmosphere.extended_T_in_K);
  auto extended_T_in_K = extended_atmosphere.extended_T_in_K;
  if (temperature_type != temperature_name) {
    Kokkos::parallel_for(
        "convert_snd2extended_atm::temperature", extended_atmos_dim, KOKKOS_LAMBDA(int i) {
          const real_t exner = i == 0 ? Kokkos::pow(p_sfc / p0, kappa)
                                      : Kokkos::pow(extended_p_in_Pa(i) / p0, kappa);
          extended_T_in_K(i) = extended_T_in_K(i) * exner;
        });
  }
  // Convert total water mixing ratio to specific humidity.
  extended_atmosphere.extended_sp_hmdty =
      read_x_profile(n_snd_var, extended_atmos_dim, rt_name, sounding_profiles);
  auto extended_sp_hmdty = extended_atmosphere.extended_sp_hmdty;
  Kokkos::parallel_for(
      "convert_snd2extended_atm::humidity", extended_atmos_dim, KOKKOS_LAMBDA(int i) {
        extended_sp_hmdty(i) = extended_sp_hmdty(i) / (extended_sp_hmdty(i) + real_t(1.));
      });
  // Read radiation scalars (currently only ozone).
  read_one_dim_file(iunit, n_rad_scalars,
                    "../input/case_data/" + sounding_trim(runtype) + "_ozone_sounding.in",
                    extended_rad_scalars);
  // New radiation-scalar arrays would be added here after ozone and in input_names.
  extended_atmosphere.extended_o3l =
      read_x_profile(1, extended_atmos_dim, ozone_name, extended_rad_scalars);
  // Free memory at the original boundary; case-owned extension arrays survive.
  extended_p_in_Pa = Real1D{};
  for (int ivar = 0; ivar < n_rad_scalars; ++ivar)
    extended_rad_scalars[ivar].values = Real1D{};
}

//=============================================================================
static void
load_extended_std_atm(int iunit,                                      // In
                      extended_atmosphere_state& extended_atmosphere) // Native module InOut
{
  //
  //  Description:
  //    Loads in the U.S. Standard atmosphere data from a file.
  //  References:
  //    McClatchey, et al., (1972) _Environmental Research Papers_,
  //    No. 411, p.94
  //
  //-----------------------------------------------------------------------------------------------

  constexpr int nCol = 5, std_atmos_dim = 50;
  const std::string atm_input_file = "../input/std_atmosphere/atmosphere.in";
  std::vector<one_dim_read_var> retVars(nCol);
  extended_atmosphere.extended_atmos_dim = std_atmos_dim;
  const int extended_atmos_dim = std_atmos_dim;
  read_one_dim_file(iunit, nCol, atm_input_file, retVars);
  // Allocate and initialize standard-atmosphere fields.
  extended_atmosphere.extended_alt =
      Real1D("load_extended_std_atm::extended_alt", extended_atmos_dim);
  extended_atmosphere.extended_T_in_K =
      Real1D("load_extended_std_atm::extended_T_in_K", extended_atmos_dim);
  extended_atmosphere.extended_sp_hmdty =
      Real1D("load_extended_std_atm::extended_sp_hmdty", extended_atmos_dim);
  extended_atmosphere.extended_p_in_mb =
      Real1D("load_extended_std_atm::extended_p_in_mb", extended_atmos_dim);
  extended_atmosphere.extended_o3l =
      Real1D("load_extended_std_atm::extended_o3l", extended_atmos_dim);
  extended_atmosphere.extended_alt =
      read_x_profile(nCol, extended_atmos_dim, z_name, retVars, atm_input_file);
  extended_atmosphere.extended_T_in_K =
      read_x_profile(nCol, extended_atmos_dim, temperature_name, retVars, atm_input_file);
  extended_atmosphere.extended_sp_hmdty =
      read_x_profile(nCol, extended_atmos_dim, sp_humidity_name, retVars, atm_input_file);
  extended_atmosphere.extended_p_in_mb =
      read_x_profile(nCol, extended_atmos_dim, press_mb_name, retVars, atm_input_file);
  extended_atmosphere.extended_o3l =
      read_x_profile(nCol, extended_atmos_dim, ozone_name, retVars, atm_input_file);
  // Deallocate memory.
  deallocate_one_dim_vars(nCol, retVars);
}

} // namespace clubb_kokkos
