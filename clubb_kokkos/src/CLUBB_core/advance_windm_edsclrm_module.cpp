// Advance mean horizontal winds and eddy scalars with Native budgets/errors.
//
// The Native Crank-Nicholson diffusion and momentum-flux treatments, perturbed
// wind solve, nudging, explicit thermodynamic copies and hole filling remain
// in source order. Physical Views are column first with LayoutLeft. Native
// matrices remain bands first until the existing solver API boundary. Scratch
// is local to each call; Native statistics owns sampled output and budgets.
// Standalone initialization rejects UV sponge damping, outside the JAX scope.
// References: Golaz et al. (2002), JAS 59, 3540–3551, Eqs. 8–9.
#include "scratch_workspace.hpp"
#include "fill_holes.hpp"
#include "mean_adv.hpp"
#include "diffusion.hpp"
#include "advance_windm_edsclrm_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "clip_explicit.hpp"
#include "advance_helper_module.hpp"
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <utility>
namespace clubb_kokkos {
namespace {
// Private named constants to avoid string comparisons.
constexpr int windm_edsclrm_um     = 1;
constexpr int windm_edsclrm_vm     = 2;
constexpr int windm_edsclrm_scalar = 3;
template <class View> std::string dump_values(View values) {
  // Native full-array diagnostic order, packed only at the host error boundary.
  auto host = to_host_space(values);
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t n = 0; n < host.size(); ++n)
    out << host.data()[n] << ' ';
  return out.str();
}
// nrhs: Number of right-hand side (explicit) vectors & Number of solution vectors.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs: Implicit contributions to um, vm, and eddy scalars  [units vary]
// rhs: Right-hand side (explicit) contributions.
// err_info: err_info struct containing err_code and err_header
// solution: Solution to the system of equations    [units vary]
template <class Array0, class Array1, class Array2>
void windm_edsclrm_solve(ScratchWorkspace& workspace,                     // InOut
                         int nzt, int ngrdcol, const grid& gr, int nrhs,  // In
                         int tridiag_solve_method,                        // In
                         bool l_implemented,                              // In
                         stats_type& stats,                               // InOut
                         bool l_need_rcond,                               // In
                         Array0 lhs, Array1 rhs, err_info_type& err_info, // InOut
                         Array2 solution)                                 // Out
    ;
// solve_type: Desc. of what is being solved for
// xm: Computed value um or vm at <t+1>    [m/s]
// invrs_dzt: The inverse spacing between momentum grid levels; centered over thermodynamic grid
// levels.
// lhs_diff: LHS diffustion terms
// lhs_ma_zt: LHS mean advection terms
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels  [m^3/kg]
// u_star_sqd: Surface friction velocity, u_star, squared      [m^2/s^2]
// rho_ds_zm: Dry, static density on momentum levels      [kg/m^3]
// wind_speed: wind speed; sqrt(u^2 + v^2)              [m/s]
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7>
void windm_edsclrm_implicit_stats(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                  int solve_type,
                                  const grid& gr,                            // In
                                  Array0 xm, Array1 invrs_dzt,               // In
                                  Array2 lhs_diff, Array3 lhs_ma_zt,         // In
                                  Array4 invrs_rho_ds_zt, Array5 u_star_sqd, // In
                                  Array6 rho_ds_zm, Array7 wind_speed,       // In
                                  bool l_imp_sfc_momentum_flux,              // In
                                  stats_type& stats)                         // InOut
    ;
// solve_type: Description of what is being solved for
// fcor: Coriolis parameter     [s^-1]
// perp_wind_m: Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
// perp_wind_g: Perpendicular component of the geostropic wind (e.g. vg)         [m/s]
// xm_forcing: Prescribed wind forcing                                          [m/s/s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// xm_tndcy: xm tendency            [m/s^2]
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void compute_uv_tndcy(ScratchWorkspace& workspace,                         // InOut
                      int nzt, int ngrdcol, int solve_type,                // In
                      Array0 fcor, Array1 perp_wind_m, Array2 perp_wind_g, // In
                      Array3 xm_forcing, bool l_implemented,               // In
                      stats_type& stats,                                   // InOut
                      Array4 xm_tndcy)                                     // Out
    ;
// dt: Model timestep                             [s]
// lhs_ma_zt: LHS mean advection terms
// lhs_diff: LHS diffustion terms
// wind_speed: wind speed; sqrt( u^2 + v^2 )              [m/s]
// u_star_sqd: Surface friction velocity, u_*, squared  [m^2/s^2]
// rho_ds_zm: Dry, static density on momentum levels     [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels [m^3/kg]
// l_implemented: True if CLUBB is being implemented and run in a host model
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
// lhs: Implicit contributions to xm (tridiagonal matrix)
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void windm_edsclrm_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, // In
                       Array0 lhs_ma_zt, Array1 lhs_diff,                        // In
                       Array2 wind_speed, Array3 u_star_sqd,                     // In
                       Array4 rho_ds_zm, Array5 invrs_rho_ds_zt,                 // In
                       bool l_implemented, bool l_imp_sfc_momentum_flux,         // In
                       Array6 lhs)                                               // Out
    ;
// solve_type: Description of what is being solved for
// dt: Model timestep                             [s]
// lhs_diff: LHS diffusion terms
// xm: Eddy-scalar variable, xm (thermo. levels)  [units vary]
// xm_tndcy: The explicit time-tendency acting on xm    [units vary]
// rho_ds_zm: Dry, static density on momentum levels     [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels [m^3/kg]
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
// xpwp_sfc: x'w' at the surface                              [units vary]
// rhs: Right-hand side (explicit) contributions.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void windm_edsclrm_rhs(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                       int solve_type,
                       real_t dt,                                     // In
                       Array0 lhs_diff, Array1 xm, Array2 xm_tndcy,   // In
                       Array3 rho_ds_zm, Array4 invrs_rho_ds_zt,      // In
                       bool l_imp_sfc_momentum_flux, Array5 xpwp_sfc, // In
                       stats_type& stats,                             // InOut
                       Array6 rhs)                                    // Out
    ;
} // namespace

//=============================================================================
// sclr_idx: Indices identifying the eddy scalars for thlm and rtm [-]
// dt: Model timestep                             [s]
// wm_zt: w wind component on thermodynamic levels    [m/s]
// Kh_zm: Eddy diffusivity for scalars on momentum levs. [m^2/s]
// clubb_params: Array of CLUBB's tunable parameters            [units vary]
// ug: u (west-to-east) geostrophic wind comp.     [m/s]
// vg: v (south-to-north) geostrophic wind comp.   [m/s]
// um_ref: Reference u wind component for nudging      [m/s]
// vm_ref: Reference v wind component for nudging      [m/s]
// wp2: w'^2 (momentum levels)                      [m^2/s^2]
// up2: u'^2 (momentum levels)                      [m^2/s^2]
// vp2: v'^2 (momentum levels)                      [m^2/s^2]
// um_forcing: u forcing                                   [m/s/s]
// vm_forcing: v forcing                                   [m/s/s]
// edsclrm_forcing: Eddy scalar large-scale forcing        [{units vary}/s]
// p_in_Pa: Air pressure (thermodynamic levels)          [Pa]
// rho_ds_zm: Dry, static density on momentum levels      [kg/m^3]
// rho_ds_zt: Dry, static density on thermodynamic levels  [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels  [m^3/kg]
// fcor: Coriolis parameter                            [s^-1]
// l_implemented: True if CLUBB is being implemented and run in a host model
// nu_vert_res_dep: Vertical resolution dependent nu values
// ts_nudge: Timescale of u/v nudging             [s]
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
// advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
// advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity when
// <u> and <v> are advanced in subroutine advance_windm_edsclrm.
// l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtm, thlm, sclrm, um and vm.
// l_uv_nudge: For wind speed nudging
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_linearize_pbl_winds: Flag (used by E3SM) to linearize PBL winds
// l_do_expldiff_rtm_thlm: Explicitly diffuse rtm and thlm with eddy scalars
// um: Mean u (west-to-east) wind component         [m/s]
// vm: Mean v (south-to-north) wind component       [m/s]
// thlm: Liquid potential temperature                 [K]
// rtm: Total water mixing ratio                     [kg/kg]
// edsclrm: Mean eddy scalar quantity             [units vary]
// upwp: <u'w'> (momentum levels)                     [m^2/s^2]
// vpwp: <v'w'> (momentum levels)                     [m^2/s^2]
// wpedsclrp: w'edsclr' (momentum levels)           [m/s {units vary}]
// um_pert: perturbed <u>    [m/s]
// vm_pert: perturbed <v>    [m/s]
// upwp_pert: perturbed <u'w'> [m^2/s^2]
// vpwp_pert: perturbed <v'w'> [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_windm_edsclrm(ScratchWorkspace& workspace,                                   // InOut
                           int nzm, int nzt, int ngrdcol, int edsclr_dim, const grid& gr, // In
                           const sclr_idx_type& sclr_idx, real_t dt,                      // In
                           Real2D wm_zt, Real2D Kh_zm, Real2D clubb_params,               // In
                           Real2D ug, Real2D vg, Real2D um_ref, Real2D vm_ref,            // In
                           Real2D wp2, Real2D up2, Real2D vp2, Real2D um_forcing,
                           Real2D vm_forcing,                                             // In
                           Real3D edsclrm_forcing, Real2D p_in_Pa,                        // In
                           Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt,    // In
                           Real1D fcor, bool l_implemented,                               // In
                           const nu_vertical_res_dep& nu_vert_res_dep, real_t ts_nudge,   // In
                           int tridiag_solve_method,                                      // In
                           bool l_predict_upwp_vpwp,                                      // In
                           bool l_upwind_xm_ma,                                           // In
                           bool l_uv_nudge,                                               // In
                           bool l_lmm_stepping,                                           // In
                           bool l_linearize_pbl_winds,                                    // In
                           bool l_do_expldiff_rtm_thlm,                                   // In
                           int fill_holes_type,                                           // In
                           int& upwp_cl_num, int& vpwp_cl_num,                            // InOut
                           stats_type& stats,                                             // InOut
                           Real2D um, Real2D vm, Real2D thlm, Real2D rtm, Real3D edsclrm, // InOut
                           Real2D upwp, Real2D vpwp, Real3D wpedsclrp,                    // InOut
                           Real2D um_pert, Real2D vm_pert, Real2D upwp_pert,
                           Real2D vpwp_pert,        // InOut
                           err_info_type& err_info) // InOut
{
  // Description:
  // Solves for both mean horizontal wind components, um and vm, and for the
  // eddy-scalars (passive scalars that don't use the high-order closure).
  // With l_do_expldiff_rtm_thlm, sclr_idx must identify distinct, valid
  // thlm/rtm slots, which check_clubb_settings_api verifies. Those slots hold
  // temporary copies for diffusion, so their edsclrm_forcing and surface
  // flux (wpedsclrp at gr%k_lb_zm, i.e. wpedsclrp_sfc) must be zero.
  // Uses the tridiagonal solver subroutine with 2 + # of scalar(s)
  // back substitutions (since the left hand side matrix is the same for all
  // input variables).
  // References:
  // Eqn. 8 & 9 on p. 3545 of
  // ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  // Method and Model Description'' Golaz, et al. (2002)
  // JAS, Vol. 59, pp. 3540--3551.
  //-----------------------------------------------------------------------
  // ------------------------ Input Variables ------------------------
  // ------------------------ Input/Output Variables ------------------------
  // Input/Output Variable for eddy-scalars
  // Output Variable for eddy-scalars
  // Variables used to track perturbed version of winds.
  // ------------------------ Local Variables ------------------------
  // Whether perturbed winds are being solved.
  auto workspace_scope = workspace.scope();
  bool l_imp_sfc_momentum_flux; // Flag for implicit momentum surface fluxes.
  bool l_need_rcond;            // Flag to request reciprocal condition number output.
  int nrhs;                     // Number of right hand side terms
  int edsclr;                   // Array index
  bool l_perturbed_wind;
  const auto sclr_idx_iiedsclr_thl = sclr_idx.iiedsclr_thl;
  const auto sclr_idx_iiedsclr_rt  = sclr_idx.iiedsclr_rt;
  const auto err_info_err_code     = err_info.err_code;
  const auto gr_k_lb_zm            = gr.k_lb_zm;
  const auto nu_vert_res_dep_nu10  = nu_vert_res_dep.nu10;
  const auto gr_k_ub_zm            = gr.k_ub_zm;

  // ------------------------ Begin Code ------------------------

  // Validate before allocating device scratch storage.
  if (edsclr_dim > 1 && l_do_expldiff_rtm_thlm) {

    // The thlm and rtm eddy scalars hold copies of thlm and rtm, which
    // advance_xm_wpxp already advances with thlm_forcing and rtm_forcing.
    // Forcing the copies here would apply those forcings a second time.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::validate_forcing", ngrdcol, KOKKOS_LAMBDA(const int column) {
          for (int k = 0; k < nzt; ++k) {
            if (edsclrm_forcing(column, k, sclr_idx_iiedsclr_thl) != zero ||
                edsclrm_forcing(column, k, sclr_idx_iiedsclr_rt) != zero)
              err_info_err_code(column) = clubb_fatal_error;
          }
        });

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "advance_windm_edsclrm: l_do_expldiff_rtm_thlm requires zero" << " "
                << " edsclrm_forcing for iiedsclr_thl and iiedsclr_rt, because" << " "
                << " thlm and rtm are already forced in advance_xm_wpxp" << "\n";
        std::cerr << message.str();
      }
      return;
    }

    // advance_xm_wpxp also already applies wpthlp_sfc and wprtp_sfc to thlm
    // and rtm. A surface flux on the copies would apply it a second time.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          if (wpedsclrp(i, gr_k_lb_zm, sclr_idx_iiedsclr_thl) != zero ||
              wpedsclrp(i, gr_k_lb_zm, sclr_idx_iiedsclr_rt) != zero) {
            err_info_err_code(i) = clubb_fatal_error;
          }
        });

    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "advance_windm_edsclrm: l_do_expldiff_rtm_thlm requires zero" << " "
                << " wpedsclrp_sfc for iiedsclr_thl and iiedsclr_rt, because" << " "
                << " wpthlp_sfc and wprtp_sfc are already applied in advance_xm_wpxp" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }
  Kokkos::View<real_t**, Kokkos::LayoutLeft> um_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::um_old", ngrdcol,
          nzt); // Saved value of mean u (west-to-east) wind component    [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vm_old =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::vm_old", ngrdcol,
          nzt); // Saved value of Mean v (south-to-north) wind component  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thlm_ed =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::thlm_ed", ngrdcol,
          nzt); // Explicit diffusion budget term for thlm
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtm_ed =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::rtm_ed", ngrdcol,
          nzt); // Explicit diffusion budget term for rtm
  Kokkos::View<real_t*, Kokkos::LayoutLeft> thlm1000 =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::thlm1000", ngrdcol);
  Kokkos::View<real_t*, Kokkos::LayoutLeft> thlm700 =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::thlm700", ngrdcol);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> edsclrm_old =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::edsclrm_old", ngrdcol, nzt,
          edsclr_dim); // Saved value of mean eddy scalar quantity   [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> um_tndcy =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::um_tndcy", ngrdcol,
          nzt); // u wind component tendency                    [m/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vm_tndcy =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::vm_tndcy", ngrdcol,
          nzt); // v wind component tendency                    [m/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_chnge =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::upwp_chnge", ngrdcol,
          nzm); // Net change of u'w' due to clipping           [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_chnge =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::vpwp_chnge", ngrdcol,
          nzm); // Net change of v'w' due to clipping           [m^2/s^2]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::lhs", 3, ngrdcol,
          nzt); // The implicit part of the tridiagonal matrix             [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::rhs", ngrdcol, nzt,
          (0 > (edsclr_dim - 2)
               ? 2
               : edsclr_dim)); // The explicit part of the tridiagonal matrix       [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> solution =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::solution", ngrdcol, nzt,
          (0 > (edsclr_dim - 2)
               ? 2
               : edsclr_dim)); // The solution to the tridiagonal matrix [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wind_speed =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::wind_speed", ngrdcol,
          nzt); // wind speed; sqrt(u^2 + v^2)              [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wind_speed_pert =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::wind_speed_pert", ngrdcol,
          nzt); // perturbed wind speed; sprt(u^2 + v^2)    [m/s]
  Kokkos::View<real_t*, Kokkos::LayoutLeft> u_star_sqd =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::u_star_sqd",
          ngrdcol); // Surface friction velocity, u_star, squared      [m^2/s^2]
  Kokkos::View<real_t*, Kokkos::LayoutLeft> u_star_sqd_pert =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::u_star_sqd_pert",
          ngrdcol); // perturbed u_star, squared                       [m^2/s^2]
  Kokkos::View<real_t*, Kokkos::LayoutLeft> nu_zero =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::nu_zero", ngrdcol);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_diff =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::lhs_diff", 3, ngrdcol,
          nzt); // LHS diffustion terms
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ma_zt =
      workspace_scope.view<Kokkos::View<real_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::lhs_ma_zt", 3, ngrdcol,
          nzt); // LHS mean advection terms
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Km_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::Km_zt", ngrdcol,
          nzt); // Eddy diffusivity of winds on thermodynamic levels [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kmh_zt =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::Kmh_zt", ngrdcol,
          nzt); // Eddy diffusivity for scalars on thermodynamic levels [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Km_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::Km_zm", ngrdcol,
          nzm); // Eddy diffusivity of winds on momentum levs. [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kmh_zm =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::Kmh_zm", ngrdcol,
          nzm); // Eddy diffusivity for scalars on momentum levels [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Km_zm_p_nu10 =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::Km_zm_p_nu10", ngrdcol,
          nzm); // Km_zm plus nu_vert_res_dep%nu10
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xpwp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "advance_windm_edsclrm::xpwp", ngrdcol,
          nzm); // x'w' for arbitrary x

  Kokkos::parallel_for(
      "advance_windm_edsclrm::loop2", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        nu_zero(i)  = zero;
      });

  Kokkos::parallel_for(
      "advance_windm_edsclrm::loop3", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = ((cell % ngrdcol));
        const int k = ((cell / ngrdcol));
        // Coefficient for momentum
        Km_zm(i, k) = Kh_zm(i, k) * clubb_params(i, ic_K10);
        // Coefficient for thermo
        Kmh_zm(i, k)       = Kh_zm(i, k) * clubb_params(i, ic_K10h);
        Km_zm_p_nu10(i, k) = Km_zm(i, k) + nu_vert_res_dep_nu10(i);
      });

  if (edsclr_dim > 1 && l_do_expldiff_rtm_thlm) {
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                          = ((cell % ngrdcol));
          const int k                          = ((cell / ngrdcol));
          edsclrm(i, k, sclr_idx_iiedsclr_thl) = thlm(i, k);
          edsclrm(i, k, sclr_idx_iiedsclr_rt)  = rtm(i, k);
        });
  }

  l_perturbed_wind = (!l_predict_upwp_vpwp) && l_linearize_pbl_winds;

  if (!l_implemented) {
    term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm, // In
                   gr.invrs_dzt, gr.invrs_dzm,                 // In
                   l_upwind_xm_ma, gr.grid_dir,                // In
                   lhs_ma_zt);                                 // Out
  } else {
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop5", std::max(0, nzt) * std::max(0, ngrdcol) * std::max(0, 3),
        KOKKOS_LAMBDA(const int cell) {
          const int j        = ((cell % 3));
          const int i        = (((cell / 3) % ngrdcol));
          const int k        = ((cell / (3 * ngrdcol)));
          lhs_ma_zt(j, i, k) = zero;
        });
  }

  if (!l_predict_upwp_vpwp) {

    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(Km_zm, Kokkos::ALL, Kokkos::ALL),
              Kokkos::subview(Km_zt, Kokkos::ALL, Kokkos::ALL), zero);

    // Calculate diffusion terms
    diffusion_zt_lhs(workspace, nzm, nzt, ngrdcol, gr, Km_zm, Km_zt, nu_vert_res_dep.nu10, // In
                     invrs_rho_ds_zt, rho_ds_zm,                                           // In
                     lhs_diff);                                                            // Out

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i  = ((cell % ngrdcol));
            const int k  = ((cell / ngrdcol));
            um_old(i, k) = um(i, k);
            vm_old(i, k) = vm(i, k);
          });
      // l_lmm_stepping
    }

    //----------------------------------------------------------------
    // Prepare tridiagonal system for horizontal winds, um and vm
    //----------------------------------------------------------------

    // Compute Coriolis, geostrophic, and other prescribed wind forcings for um.
    compute_uv_tndcy(workspace, nzt, ngrdcol, windm_edsclrm_um, // In
                     fcor, vm, vg,                              // In
                     um_forcing, l_implemented,                 // In
                     stats,                                     // InOut
                     um_tndcy);                                 // Out

    // Compute Coriolis, geostrophic, and other prescribed wind forcings for vm.
    compute_uv_tndcy(workspace, nzt, ngrdcol, windm_edsclrm_vm, // In
                     fcor, um, ug,                              // In
                     vm_forcing, l_implemented,                 // In
                     stats,                                     // InOut
                     vm_tndcy);                                 // Out

    // Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied through
    // an implicit method, such that:
    //    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
    l_imp_sfc_momentum_flux = true;

    // Compute wind speed (use threshold "eps" to prevent divide-by-zero
    // error).
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop7", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i      = ((cell % ngrdcol));
          const int k      = ((cell / ngrdcol));
          wind_speed(i, k) = ((Kokkos::sqrt((um(i, k) * um(i, k)) + (vm(i, k) * vm(i, k)))) > (eps)
                                  ? (Kokkos::sqrt((um(i, k) * um(i, k)) + (vm(i, k) * vm(i, k))))
                                  : (eps));
        });

    // Compute u_star_sqd according to the definition of u_star.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop8", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i   = cell;
          u_star_sqd(i) = Kokkos::sqrt((upwp(i, gr_k_lb_zm) * upwp(i, gr_k_lb_zm)) +
                                       (vpwp(i, gr_k_lb_zm) * vpwp(i, gr_k_lb_zm)));
        });

    // Compute the explicit portion of the um equation.
    // Build the right-hand side vector.
    windm_edsclrm_rhs(
        workspace, nzm, nzt, ngrdcol, gr, windm_edsclrm_um, dt,                    // In
        lhs_diff, um, um_tndcy,                                                    // In
        rho_ds_zm, invrs_rho_ds_zt,                                                // In
        l_imp_sfc_momentum_flux, Kokkos::subview(upwp, Kokkos::ALL, (gr.k_lb_zm)), // In
        stats,                                                                     // InOut
        Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, ((windm_edsclrm_um - 1)))); // Out

    // Compute the explicit portion of the vm equation.
    // Build the right-hand side vector.
    windm_edsclrm_rhs(
        workspace, nzm, nzt, ngrdcol, gr, windm_edsclrm_vm, dt,                    // In
        lhs_diff, vm, vm_tndcy,                                                    // In
        rho_ds_zm, invrs_rho_ds_zt,                                                // In
        l_imp_sfc_momentum_flux, Kokkos::subview(vpwp, Kokkos::ALL, (gr.k_lb_zm)), // In
        stats,                                                                     // InOut
        Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, ((windm_edsclrm_vm - 1)))); // Out

    // Store momentum flux (explicit component)
    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, um,      // In
              xpwp);                 // Out

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop9", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol) + 1);
          upwp(i, k)  = -one_half * xpwp(i, k);
        });

    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, vm,      // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop10", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol) + 1);
          vpwp(i, k)  = -one_half * xpwp(i, k);
        });

    // A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
    // means that x'w' at the top model level is 0,
    // since x'w' = - K_zm * d(xm)/dz.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop11", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i         = cell;
          upwp(i, gr_k_ub_zm) = zero;
          vpwp(i, gr_k_ub_zm) = zero;
        });

    // Compute the implicit portion of the um and vm equations.
    // Build the left-hand side matrix.
    windm_edsclrm_lhs(nzm, nzt, ngrdcol, gr, dt,              // In
                      lhs_ma_zt, lhs_diff,                    // In
                      wind_speed, u_star_sqd,                 // In
                      rho_ds_zm, invrs_rho_ds_zt,             // In
                      l_implemented, l_imp_sfc_momentum_flux, // In
                      lhs);                                   // Out

    // Decompose and back substitute for um and vm
    nrhs         = 2;
    l_need_rcond = stats.l_sample && var_on_stats_list(stats, "windm_matrix_condt_num");
    windm_edsclrm_solve(workspace, nzt, ngrdcol, gr, nrhs, // In
                        tridiag_solve_method,              // In
                        l_implemented,                     // In
                        stats,                             // InOut
                        l_need_rcond,                      // In
                        lhs, rhs, err_info,                // InOut
                        solution);                         // Out

    // Check for singular matrices and bad LAPACK arguments
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "Fatal error solving for um/vm" << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }

    //----------------------------------------------------------------
    // Update zonal (west-to-east) component of mean wind, um
    //----------------------------------------------------------------
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop12", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol));
          um(i, k)    = solution(i, k, (windm_edsclrm_um - 1));
        });

    //----------------------------------------------------------------
    // Update meridional (south-to-north) component of mean wind, vm
    //----------------------------------------------------------------
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop13", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol));
          vm(i, k)    = solution(i, k, (windm_edsclrm_vm - 1));
        });

    if (stats.l_sample) {

      // Implicit contributions to um and vm
      windm_edsclrm_implicit_stats(workspace, nzm, nzt, ngrdcol, windm_edsclrm_um, gr, // In
                                   um, gr.invrs_dzt,                                   // In
                                   lhs_diff, lhs_ma_zt,                                // In
                                   invrs_rho_ds_zt, u_star_sqd,                        // In
                                   rho_ds_zm, wind_speed,                              // In
                                   l_imp_sfc_momentum_flux,                            // In
                                   stats);                                             // InOut

      windm_edsclrm_implicit_stats(workspace, nzm, nzt, ngrdcol, windm_edsclrm_vm, gr, // In
                                   vm, gr.invrs_dzt,                                   // In
                                   lhs_diff, lhs_ma_zt,                                // In
                                   invrs_rho_ds_zt, u_star_sqd,                        // In
                                   rho_ds_zm, wind_speed,                              // In
                                   l_imp_sfc_momentum_flux,                            // In
                                   stats);                                             // InOut
    }

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop14", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            um(i, k)    = one_half * (um_old(i, k) + um(i, k));
            vm(i, k)    = one_half * (vm_old(i, k) + vm(i, k));
          });
      // l_lmm_stepping ) then
    }

    // Native uv_sponge_damp_settings%l_sponge_damping:
    // Sponge damping is excluded with the JAX scope.
    // Standalone initialization rejects sponge-enabled cases before allocation.

    // Second part of momentum (implicit component)

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used.
    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, um,      // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop15", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol) + 1);
          upwp(i, k)  = upwp(i, k) - one_half * xpwp(i, k);
        });

    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, vm,      // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop16", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol) + 1);
          vpwp(i, k)  = vpwp(i, k) - one_half * xpwp(i, k);
        });

    // Adjust um and vm if nudging is turned on.
    if (l_uv_nudge) {

      // Reflect nudging in budget
      if (stats.l_sample) {
        {
          Real2D stats_values_17 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_17", um.extent(0),
              um.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", um.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % um.extent(0)));
                const int k           = ((cell / um.extent(0)));
                stats_values_17(i, k) = um(i, k) / dt;
              });
          stats_begin_budget("um_ndg", stats_values_17, // In
                             stats);                    // InOut
        }
        {
          Real2D stats_values_18 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_18", vm.extent(0),
              vm.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vm.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vm.extent(0)));
                const int k           = ((cell / vm.extent(0)));
                stats_values_18(i, k) = vm(i, k) / dt;
              });
          stats_begin_budget("vm_ndg", stats_values_18, // In
                             stats);                    // InOut
        }
      }

      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop19", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            um(i, k)    = um(i, k) - ((um(i, k) - um_ref(i, k)) * (dt / ts_nudge));
            vm(i, k)    = vm(i, k) - ((vm(i, k) - vm_ref(i, k)) * (dt / ts_nudge));
          });

      if (stats.l_sample) {
        {
          Real2D stats_values_20 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_20", um.extent(0),
              um.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", um.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % um.extent(0)));
                const int k           = ((cell / um.extent(0)));
                stats_values_20(i, k) = um(i, k) / dt;
              });
          stats_finalize_budget("um_ndg", stats_values_20, // In
                                stats);                    // InOut
        }
        {
          Real2D stats_values_21 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_21", vm.extent(0),
              vm.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vm.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vm.extent(0)));
                const int k           = ((cell / vm.extent(0)));
                stats_values_21(i, k) = vm(i, k) / dt;
              });
          stats_finalize_budget("vm_ndg", stats_values_21, // In
                                stats);                    // InOut
        }
      }

      // l_uv_nudge
    }

    if (stats.l_sample) {
      stats_update("um_ref", um_ref, // In
                   stats);           // InOut
      stats_update("vm_ref", vm_ref, // In
                   stats);           // InOut
    }

    // Clipping for u'w'
    //
    // Clipping u'w' at each vertical level, based on the
    // correlation of u and w at each vertical level, such that:
    // corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
    // -1 <= corr_(u,w) <= 1.
    //
    // Since u'^2, w'^2, and u'w' are each advanced in different
    // subroutines from each other in advance_clubb_core, clipping for u'w'
    // has to be done three times during each timestep (once after each
    // variable has been updated).
    // This clip can be the first, middle, or last budget contribution,
    // depending on the relative advancement order in this timestep.
    if (stats.l_sample && l_predict_upwp_vpwp) {
      if (upwp_cl_num == 0) {
        {
          Real2D stats_values_22 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_22", upwp.extent(0),
              upwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % upwp.extent(0)));
                const int k           = ((cell / upwp.extent(0)));
                stats_values_22(i, k) = upwp(i, k) / dt;
              });
          stats_begin_budget("upwp_cl", stats_values_22, // In
                             stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_23 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_23", upwp.extent(0),
              upwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % upwp.extent(0)));
                const int k           = ((cell / upwp.extent(0)));
                stats_values_23(i, k) = -upwp(i, k) / dt;
              });
          stats_update_budget("upwp_cl", stats_values_23, // In
                              stats);                     // InOut
        }
      }
    }
    upwp_cl_num = (upwp_cl_num + 1);
    clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, // In
               upwp,                              // InOut
               upwp_chnge);                       // Out
    if (stats.l_sample && l_predict_upwp_vpwp) {
      if (upwp_cl_num == upwp_cl_max) {
        {
          Real2D stats_values_24 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_24", upwp.extent(0),
              upwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % upwp.extent(0)));
                const int k           = ((cell / upwp.extent(0)));
                stats_values_24(i, k) = upwp(i, k) / dt;
              });
          stats_finalize_budget("upwp_cl", stats_values_24, // In
                                stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_25 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_25", upwp.extent(0),
              upwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", upwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % upwp.extent(0)));
                const int k           = ((cell / upwp.extent(0)));
                stats_values_25(i, k) = upwp(i, k) / dt;
              });
          stats_update_budget("upwp_cl", stats_values_25, // In
                              stats);                     // InOut
        }
      }
    }

    // Clipping for v'w'
    //
    // Clipping v'w' at each vertical level, based on the
    // correlation of v and w at each vertical level, such that:
    // corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
    // -1 <= corr_(v,w) <= 1.
    //
    // Since v'^2, w'^2, and v'w' are each advanced in different
    // subroutines from each other in advance_clubb_core, clipping for v'w'
    // has to be done three times during each timestep (once after each
    // variable has been updated).
    // This clip can be the first, middle, or last budget contribution,
    // depending on the relative advancement order in this timestep.
    if (stats.l_sample && l_predict_upwp_vpwp) {
      if (vpwp_cl_num == 0) {
        {
          Real2D stats_values_26 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_26", vpwp.extent(0),
              vpwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vpwp.extent(0)));
                const int k           = ((cell / vpwp.extent(0)));
                stats_values_26(i, k) = vpwp(i, k) / dt;
              });
          stats_begin_budget("vpwp_cl", stats_values_26, // In
                             stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_27 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_27", vpwp.extent(0),
              vpwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vpwp.extent(0)));
                const int k           = ((cell / vpwp.extent(0)));
                stats_values_27(i, k) = -vpwp(i, k) / dt;
              });
          stats_update_budget("vpwp_cl", stats_values_27, // In
                              stats);                     // InOut
        }
      }
    }
    vpwp_cl_num = (vpwp_cl_num + 1);
    clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, // In
               vpwp,                              // InOut
               vpwp_chnge);                       // Out
    if (stats.l_sample && l_predict_upwp_vpwp) {
      if (vpwp_cl_num == vpwp_cl_max) {
        {
          Real2D stats_values_28 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_28", vpwp.extent(0),
              vpwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vpwp.extent(0)));
                const int k           = ((cell / vpwp.extent(0)));
                stats_values_28(i, k) = vpwp(i, k) / dt;
              });
          stats_finalize_budget("vpwp_cl", stats_values_28, // In
                                stats);                     // InOut
        }
      } else {
        {
          Real2D stats_values_29 = workspace_scope.view<Real2D>(
              Kokkos::WithoutInitializing, "advance_windm_edsclrm::stats_values_29", vpwp.extent(0),
              vpwp.extent(1));
          Kokkos::parallel_for(
              "advance_windm_edsclrm::stats_value", vpwp.size(), KOKKOS_LAMBDA(int cell) {
                const int i           = ((cell % vpwp.extent(0)));
                const int k           = ((cell / vpwp.extent(0)));
                stats_values_29(i, k) = vpwp(i, k) / dt;
              });
          stats_update_budget("vpwp_cl", stats_values_29, // In
                              stats);                     // InOut
        }
      }
    }

    // .not. l_predict_upwp_vpwp
  }

  if (l_perturbed_wind) {

    //----------------------------------------------------------------
    // Prepare tridiagonal system for horizontal winds, um and vm
    //----------------------------------------------------------------

    // Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied through
    // an implicit method, such that:
    //    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
    l_imp_sfc_momentum_flux = true;

    // Compute wind speed (use threshold "eps" to prevent divide-by-zero
    // error).
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop30", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i           = ((cell % ngrdcol));
          const int k           = ((cell / ngrdcol));
          wind_speed_pert(i, k) = ((Kokkos::sqrt(((um_pert(i, k)) * (um_pert(i, k))) +
                                                 ((vm_pert(i, k)) * (vm_pert(i, k))))) > (eps)
                                       ? (Kokkos::sqrt(((um_pert(i, k)) * (um_pert(i, k))) +
                                                       ((vm_pert(i, k)) * (vm_pert(i, k)))))
                                       : (eps));
        });

    // Compute u_star_sqd according to the definition of u_star.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop31", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i        = cell;
          u_star_sqd_pert(i) = Kokkos::sqrt((upwp_pert(i, gr_k_lb_zm) * upwp_pert(i, gr_k_lb_zm)) +
                                            (vpwp_pert(i, gr_k_lb_zm) * vpwp_pert(i, gr_k_lb_zm)));
        });

    // Compute the explicit portion of the um equation.
    // Build the right-hand side vector.
    windm_edsclrm_rhs(
        workspace, nzm, nzt, ngrdcol, gr, windm_edsclrm_um, dt,                         // In
        lhs_diff, um_pert, um_tndcy,                                                    // In
        rho_ds_zm, invrs_rho_ds_zt,                                                     // In
        l_imp_sfc_momentum_flux, Kokkos::subview(upwp_pert, Kokkos::ALL, (gr.k_lb_zm)), // In
        stats,                                                                          // InOut
        Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, ((windm_edsclrm_um - 1))));      // Out

    // Compute the explicit portion of the vm equation.
    // Build the right-hand side vector.
    windm_edsclrm_rhs(
        workspace, nzm, nzt, ngrdcol, gr, windm_edsclrm_vm, dt,                         // In
        lhs_diff, vm_pert, vm_tndcy,                                                    // In
        rho_ds_zm, invrs_rho_ds_zt,                                                     // In
        l_imp_sfc_momentum_flux, Kokkos::subview(vpwp_pert, Kokkos::ALL, (gr.k_lb_zm)), // In
        stats,                                                                          // InOut
        Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, ((windm_edsclrm_vm - 1))));      // Out

    // Store momentum flux (explicit component)

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used.
    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, um_pert, // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop32", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i     = ((cell % ngrdcol));
          const int k     = ((cell / ngrdcol) + 1);
          upwp_pert(i, k) = -one_half * xpwp(i, k);
        });

    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, vm_pert, // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop33", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i     = ((cell % ngrdcol));
          const int k     = ((cell / ngrdcol) + 1);
          vpwp_pert(i, k) = -one_half * xpwp(i, k);
        });

    // A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
    // means that x'w' at the top model level is 0,
    // since x'w' = - K_zm * d(xm)/dz.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop34", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i              = cell;
          upwp_pert(i, gr_k_ub_zm) = zero;
          vpwp_pert(i, gr_k_ub_zm) = zero;
        });

    // Compute the implicit portion of the um and vm equations.
    // Build the left-hand side matrix.
    windm_edsclrm_lhs(nzm, nzt, ngrdcol, gr, dt,              // In
                      lhs_ma_zt, lhs_diff,                    // In
                      wind_speed_pert, u_star_sqd_pert,       // In
                      rho_ds_zm, invrs_rho_ds_zt,             // In
                      l_implemented, l_imp_sfc_momentum_flux, // In
                      lhs);                                   // Out

    // Decompose and back substitute for um and vm
    nrhs         = 2;
    l_need_rcond = stats.l_sample && var_on_stats_list(stats, "windm_matrix_condt_num");
    windm_edsclrm_solve(workspace, nzt, ngrdcol, gr, nrhs, // In
                        tridiag_solve_method,              // In
                        l_implemented,                     // In
                        stats,                             // InOut
                        l_need_rcond,                      // In
                        lhs, rhs, err_info,                // InOut
                        solution);                         // Out

    // Check for singular matrices and bad LAPACK arguments
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "Fatal error solving for um_pert/vm_pert" << "\n";
          std::cerr << message.str();
        }
        return;
      }
    }

    //----------------------------------------------------------------
    // Update zonal (west-to-east) component of mean wind, um
    //----------------------------------------------------------------
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop35", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i   = ((cell % ngrdcol));
          const int k   = ((cell / ngrdcol));
          um_pert(i, k) = solution(i, k, (windm_edsclrm_um - 1));
        });

    //----------------------------------------------------------------
    // Update meridional (south-to-north) component of mean wind, vm
    //----------------------------------------------------------------
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop36", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i   = ((cell % ngrdcol));
          const int k   = ((cell / ngrdcol));
          vm_pert(i, k) = solution(i, k, (windm_edsclrm_vm - 1));
        });

    // Second part of momentum (implicit component)

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used.
    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, um_pert, // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop37", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i     = ((cell % ngrdcol));
          const int k     = ((cell / ngrdcol) + 1);
          upwp_pert(i, k) = upwp_pert(i, k) - one_half * xpwp(i, k);
        });

    calc_xpwp(nzm, nzt, ngrdcol, gr, // In
              Km_zm_p_nu10, vm_pert, // In
              xpwp);                 // Out

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop38", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i     = ((cell % ngrdcol));
          const int k     = ((cell / ngrdcol) + 1);
          vpwp_pert(i, k) = vpwp_pert(i, k) - one_half * xpwp(i, k);
        });

    // Clipping for u'w'
    //
    // Clipping u'w' at each vertical level, based on the
    // correlation of u and w at each vertical level, such that:
    // corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
    // -1 <= corr_(u,w) <= 1.
    //
    // Since u'^2, w'^2, and u'w' are each advanced in different
    // subroutines from each other in advance_clubb_core, clipping for u'w'
    // has to be done three times during each timestep (once after each
    // variable has been updated).
    // This is the third instance of u'w' clipping.
    clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, // In
               upwp_pert,                         // InOut
               upwp_chnge);                       // Out

    // Clipping for v'w'
    //
    // Clipping v'w' at each vertical level, based on the
    // correlation of v and w at each vertical level, such that:
    // corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
    // -1 <= corr_(v,w) <= 1.
    //
    // Since v'^2, w'^2, and v'w' are each advanced in different
    // subroutines from each other in advance_clubb_core, clipping for v'w'
    // has to be done three times during each timestep (once after each
    // variable has been updated).
    // This is the third instance of v'w' clipping.
    clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, // In
               vpwp_pert,                         // InOut
               vpwp_chnge);                       // Out

    // l_perturbed_wind
  }

  //----------------------------------------------------------------
  // Prepare tridiagonal system for eddy-scalars
  //----------------------------------------------------------------

  if (edsclr_dim > 0) {

    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(Kmh_zm, Kokkos::ALL, Kokkos::ALL),
              Kokkos::subview(Kmh_zt, Kokkos::ALL, Kokkos::ALL), zero);

    // Calculate diffusion terms
    diffusion_zt_lhs(workspace, nzm, nzt, ngrdcol, gr, Kmh_zm, Kmh_zt, nu_zero, // In
                     invrs_rho_ds_zt, rho_ds_zm,                                // In
                     lhs_diff);                                                 // Out

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop39",
          std::max(0, edsclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i               = ((cell % ngrdcol));
            const int k               = (((cell / ngrdcol) % nzt));
            const int edsclr          = ((cell / (ngrdcol * nzt)));
            edsclrm_old(i, k, edsclr) = edsclrm(i, k, edsclr);
          });
      // l_lmm_stepping
    }

    // Eddy-scalar surface fluxes, x'w'|_sfc, are applied through an explicit
    // method.
    l_imp_sfc_momentum_flux = false;

    // Compute the explicit portion of eddy scalar equation.
    // Build the right-hand side vector.
    // Because of statistics, we have to use a DO rather than a FORALL here
    // -dschanen 7 Oct 2008
    for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
      windm_edsclrm_rhs(workspace, nzm, nzt, ngrdcol, gr,                                     // In
                        windm_edsclrm_scalar, dt,                                             // In
                        lhs_diff, Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr), // In
                        Kokkos::subview(edsclrm_forcing, Kokkos::ALL, Kokkos::ALL, edsclr),   // In
                        rho_ds_zm, invrs_rho_ds_zt,                                           // In
                        l_imp_sfc_momentum_flux,                                              // In
                        Kokkos::subview(wpedsclrp, Kokkos::ALL, (gr.k_lb_zm), edsclr),        // In
                        stats,                                                   // InOut
                        Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, edsclr)); // Out
    }

    // Store momentum flux (explicit component)

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used
    for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {

      calc_xpwp(nzm, nzt, ngrdcol, gr,                                                    // In
                Km_zm_p_nu10, Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr), // In
                xpwp);                                                                    // Out

      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop40", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i             = ((cell % ngrdcol));
            const int k             = ((cell / ngrdcol) + 1);
            wpedsclrp(i, k, edsclr) = -one_half * xpwp(i, k);
          });
    }

    // A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
    // means that x'w' at the top model level is 0,
    // since x'w' = - K_zm * d(xm)/dz.
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop41", std::max(0, edsclr_dim) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i                      = ((cell % ngrdcol));
          const int edsclr                 = ((cell / ngrdcol));
          wpedsclrp(i, gr_k_ub_zm, edsclr) = zero;
        });

    // Compute the implicit portion of the xm (eddy-scalar) equations.
    // Build the left-hand side matrix.
    windm_edsclrm_lhs(nzm, nzt, ngrdcol, gr, dt,              // In
                      lhs_ma_zt, lhs_diff,                    // In
                      wind_speed, u_star_sqd,                 // In
                      rho_ds_zm, invrs_rho_ds_zt,             // In
                      l_implemented, l_imp_sfc_momentum_flux, // In
                      lhs);                                   // Out

    // Decompose and back substitute for all eddy-scalar variables
    l_need_rcond = false;
    windm_edsclrm_solve(workspace, nzt, ngrdcol, gr, edsclr_dim, // In
                        tridiag_solve_method,                    // In
                        l_implemented,                           // In
                        stats,                                   // InOut
                        l_need_rcond,                            // In
                        lhs, rhs, err_info,                      // InOut
                        solution);                               // Out

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << err_info.err_header_global << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "Fatal error solving for edsclrm" << "\n";
          std::cerr << message.str();
        }
      }
    }

    //----------------------------------------------------------------
    // Update Eddy-diff. Passive Scalars
    //----------------------------------------------------------------
    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop42",
        std::max(0, edsclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i           = ((cell % ngrdcol));
          const int k           = (((cell / ngrdcol) % nzt));
          const int edsclr      = ((cell / (ngrdcol * nzt)));
          edsclrm(i, k, edsclr) = solution(i, k, edsclr);
        });

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop43",
          std::max(0, edsclr_dim) * std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i           = ((cell % ngrdcol));
            const int k           = (((cell / ngrdcol) % nzt));
            const int edsclr      = ((cell / (ngrdcol * nzt)));
            edsclrm(i, k, edsclr) = one_half * (edsclrm_old(i, k, edsclr) + edsclrm(i, k, edsclr));
          });
      // l_lmm_stepping
    }

    // Second part of momentum (implicit component)

    // Solve for x'w' at all intermediate model levels.
    // A Crank-Nicholson timestep is used.
    for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {

      calc_xpwp(nzm, nzt, ngrdcol, gr,                                              // In
                Kmh_zm, Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr), // In
                xpwp);                                                              // Out

      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop44", std::max(0, (nzm - 2)) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i             = ((cell % ngrdcol));
            const int k             = ((cell / ngrdcol) + 1);
            wpedsclrp(i, k, edsclr) = -one_half * xpwp(i, k);
          });
    }

    // Note that the w'edsclr' terms are not clipped, since we don't compute
    // the variance of edsclr anywhere. -dschanen 7 Oct 2008
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {

      // err_header already printed earlier
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Error in advance_windm_edsclrm" << "\n";
        std::cerr << message.str();
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "intent(in)" << "\n";
        std::cerr << message.str();
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "dt = " << " " << dt << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wm_zt = " << " " << dump_values(wm_zt) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Km_zm = " << " " << dump_values(Km_zm) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "ug = " << " " << dump_values(ug) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vg = " << " " << dump_values(vg) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "um_ref = " << " " << dump_values(um_ref) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vm_ref = " << " " << dump_values(vm_ref) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wp2 = " << " " << dump_values(wp2) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "up2 = " << " " << dump_values(up2) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vp2 = " << " " << dump_values(vp2) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "um_forcing = " << " " << dump_values(um_forcing) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vm_forcing = " << " " << dump_values(vm_forcing) << "\n";
        std::cerr << message.str();
      }
      for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "edsclrm_forcing # = " << " " << ((edsclr + 1)) << " "
                  << dump_values(edsclrm_forcing) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "fcor = " << " " << dump_values(fcor) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "l_implemented = " << " " << l_implemented << "\n";
        std::cerr << message.str();
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "intent(inout)" << "\n";
        std::cerr << message.str();
      }

      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "um (pre-solve) = " << " " << dump_values(um_old) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "um = " << " " << dump_values(um) << "\n";
        std::cerr << message.str();
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "vm (pre-solve) = " << " " << dump_values(vm_old) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vm = " << " " << dump_values(vm) << "\n";
        std::cerr << message.str();
      }
      for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
        if (l_lmm_stepping) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "edsclrm (pre-solve) # " << " " << ((edsclr + 1)) << " " << "=" << " "
                    << dump_values(Kokkos::subview(edsclrm_old, Kokkos::ALL, Kokkos::ALL, edsclr))
                    << "\n";
            std::cerr << message.str();
          }
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "edsclrm # " << " " << ((edsclr + 1)) << " " << "=" << " "
                  << dump_values(Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr))
                  << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "upwp = " << " " << dump_values(upwp) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vpwp = " << " " << dump_values(vpwp) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wpedsclrp = " << " " << dump_values(wpedsclrp) << "\n";
        std::cerr << message.str();
      }

      return;
    }
  }

  if (edsclr_dim > 1 && l_do_expldiff_rtm_thlm) {

    pvertinterp(nzt, ngrdcol, gr,               // In
                p_in_Pa, real_t(70000.0), thlm, // In
                thlm700);                       // Out

    pvertinterp(nzt, ngrdcol, gr,                // In
                p_in_Pa, real_t(100000.0), thlm, // In
                thlm1000);                       // Out

    if (stats.l_sample) {

      // thlm_ed and rtm_ed are budget terms intended to track the effect of explicit diffusion
      Kokkos::parallel_for(
          "advance_windm_edsclrm::loop45", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            if (thlm700(i) - thlm1000(i) < real_t(20.0)) {
              thlm_ed(i, k) = (edsclrm(i, k, sclr_idx_iiedsclr_thl) - thlm(i, k)) / dt;
              rtm_ed(i, k)  = (edsclrm(i, k, sclr_idx_iiedsclr_rt) - rtm(i, k)) / dt;
            } else {
              thlm_ed(i, k) = zero;
              rtm_ed(i, k)  = zero;
            }
          });
    }

    Kokkos::parallel_for(
        "advance_windm_edsclrm::loop46", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = ((cell % ngrdcol));
          const int k = ((cell / ngrdcol));
          if (thlm700(i) - thlm1000(i) < real_t(20.0)) {
            thlm(i, k) = edsclrm(i, k, sclr_idx_iiedsclr_thl);
            rtm(i, k)  = edsclrm(i, k, sclr_idx_iiedsclr_rt);
          }
        });
  }

  if (stats.l_sample) {
    if (edsclr_dim < 2 || !l_do_expldiff_rtm_thlm) {
      Kokkos::parallel_for(
          "advance_windm_edsclrm::array47", (thlm_ed.extent(0) * thlm_ed.extent(1)),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0         = ((cell % thlm_ed.extent(0)));
            const int slice_1         = ((cell / thlm_ed.extent(0)));
            thlm_ed(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "advance_windm_edsclrm::array48", (rtm_ed.extent(0) * rtm_ed.extent(1)),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0        = ((cell % rtm_ed.extent(0)));
            const int slice_1        = ((cell / rtm_ed.extent(0)));
            rtm_ed(slice_0, slice_1) = zero;
          });
    }
    stats_update("thlm_ed", thlm_ed, // In
                 stats);             // InOut
    stats_update("rtm_ed", rtm_ed,   // In
                 stats);             // InOut
  }

  if (edsclr_dim > 0 && fill_holes_type != 0) {
    for (edsclr = 0; edsclr < edsclr_dim; edsclr += 1) {
      fill_holes_vertical_api(workspace, nzt, ngrdcol, zero_threshold,                     // In
                              gr.k_lb_zt, gr.k_ub_zt,                                      // In
                              gr.dzt, rho_ds_zt, gr.grid_dir_indx,                         // In
                              fill_holes_type,                                             // In
                              Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr)); // InOut
    }
  }

  return;
}

namespace {
//=============================================================================
// nrhs: Number of right-hand side (explicit) vectors & Number of solution vectors.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// l_implemented: True if CLUBB is being implemented and run in a host model
// lhs: Implicit contributions to um, vm, and eddy scalars  [units vary]
// rhs: Right-hand side (explicit) contributions.
// err_info: err_info struct containing err_code and err_header
// solution: Solution to the system of equations    [units vary]
template <class Array0, class Array1, class Array2>
void windm_edsclrm_solve(ScratchWorkspace& workspace,                     // InOut
                         int nzt, int ngrdcol, const grid& gr, int nrhs,  // In
                         int tridiag_solve_method,                        // In
                         bool l_implemented,                              // In
                         stats_type& stats,                               // InOut
                         bool l_need_rcond,                               // In
                         Array0 lhs, Array1 rhs, err_info_type& err_info, // InOut
                         Array2 solution)                                 // Out
{
  // Description:
  // Solves the horizontal wind or eddy-scalar time-tendency equation, and
  // diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm
  // is used in solving the turbulent advection term and in diagnosing the
  // turbulent flux.
  //
  // The rate of change of an eddy-scalar quantity, xm, is:
  //
  // d(xm)/dt = - w * d(xm)/dz - (1/rho_ds) * d( rho_ds * x'w' )/dz
  //            + xm_forcings.
  //
  //
  // The Turbulent Advection Term
  // ----------------------------
  //
  // The above equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * x'w' )/dz;
  //
  // where the momentum flux, x'w', is closed using a down gradient approach:
  //
  // x'w' = - K_zm * d(xm)/dz.
  //
  // The turbulent advection term becomes:
  //
  // + (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz;
  //
  // which is the same as a standard eddy-diffusion term (if "rho_ds * K_zm" in
  // the term above is substituted for "K_zm" in a standard eddy-diffusion
  // term, and if the standard eddy-diffusion term is multiplied by
  // "1/rho_ds").  Thus, the turbulent advection term is treated and solved in
  // the same way that a standard eddy-diffusion term would be solved.  The
  // term is discretized as follows:
  //
  // The values of xm are found on the thermodynamic levels, while the values
  // of K_zm are found on the momentum levels.  Additionally, the values of
  // rho_ds_zm are found on the momentum levels, and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  The
  // derivatives (d/dz) of xm are taken over the intermediate momentum levels.
  // At the intermediate momentum levels, d(xm)/dz is multiplied by K_zm and by
  // rho_ds_zm.  Then, the derivative of the whole mathematical expression is
  // taken over the central thermodynamic level, where it is multiplied by
  // invrs_rho_ds_zt, which yields the desired result.
  //
  // ---xm(kp1)----------------------------------------------------- t(k+1)
  //
  // ===========d(xm)/dz===K_zm(k+1)===rho_ds_zm(k+1)=============== m(k+1)
  //
  // ---xm(k)---invrs_rho_ds_zt---d[rho_ds_zm*K_zm*d(xm)/dz]/dz----- t(k)
  //
  // ===========d(xm)/dz===K_zm(k)=====rho_ds_zm(k)================= m(k)
  //
  // ---xm(km1)----------------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is used
  // for momentum levels.
  //
  // invrs_dzt(k)   = 1 / ( zm(k+1) - zm(k) )
  // invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
  // invrs_dzm(k)   = 1 / ( zt(k) - zt(k-1) )
  //
  // The vertically discretized form of the turbulent advection term (treated
  // as an eddy diffusion term) is written out as:
  //
  // + invrs_rho_ds_zt(k)
  //   * invrs_dzt(k)
  //     * [   rho_ds_zm(k+1) * K_zm(k+1)
  //           * invrs_dzm(k+1) * ( xm(k+1) - xm(k) )
  //         - rho_ds_zm(k) * K_zm(k)
  //           * invrs_dzm(k) * ( xm(k) - xm(k-1) ) ].
  //
  // For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is
  // used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
  // eddy-diffusion term.  The discretized implicit portion of the term is
  // written out as:
  //
  // + (1/2) * invrs_rho_ds_zt(k)
  //   * invrs_dzt(k)
  //     * [   rho_ds_zm(k+1) * K_zm(k+1)
  //           * invrs_dzm(k+1) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
  //         - rho_ds_zm(k) * K_zm(k)
  //           * invrs_dzm(k) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].
  //
  // Note:  When the implicit term is brought over to the left-hand side,
  //        the sign is reversed and the leading "+" in front of the term
  //        is changed to a "-".
  //
  // The discretized explicit portion of the term is written out as:
  //
  // + (1/2) * invrs_rho_ds_zt(k)
  //   * invrs_dzt(k)
  //     * [   rho_ds_zm(k+1) * K_zm(k+1)
  //           * invrs_dzm(k+1) * ( xm(k+1,<t>) - xm(k,<t>) )
  //         - rho_ds_zm(k) * K_zm(k)
  //           * invrs_dzm(k) * ( xm(k,<t>) - xm(k-1,<t>) ) ].
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(xm)/dt equation.
  //
  //
  // Boundary Conditions:
  //
  // An eddy-scalar quantity is not allowed to flux out the upper boundary.
  // Thus, a zero-flux boundary condition is used for the upper boundary in the
  // eddy-diffusion equation.
  //
  // The lower boundary condition is much more complicated.  It is neither a
  // zero-flux nor a fixed-point boundary condition.  Rather, it is a
  // fixed-flux boundary condition.  This term is a turbulent advection term,
  // but with the eddy-scalars, the only value of x'w' relevant in solving the
  // d(xm)/dt equation is the value of x'w' at the surface (the first momentum
  // level), which is written as x'w'|_sfc.
  //
  // 1) x'w' surface flux; generalized explicit form
  //
  //    The x'w' surface flux is applied to the d(xm)/dt equation through the
  //    turbulent advection term, which is:
  //
  //    - (1/rho_ds) * d( rho_ds * x'w' )/dz.
  //
  //    At most vertical levels, a substitution can be made for x'w', such
  //    that:
  //
  //    x'w' = - K_zm * d(xm)/dz.
  //
  //    However, the same substitution cannot be made at the surface (momentum
  //    level 1), as x'w'|_sfc is a surface flux that is explicitly computed
  //    elsewhere in the model code.
  //
  //    The lower boundary condition, which in this case needs to be applied to
  //    the d(xm)/dt equation at level 1, is discretized as follows:
  //
  //    --xm(2)------------------------------------------------------- t(2)
  //
  //    ========[x'w'(2) = -K_zm(2)*d(xm)/dz]===rho_ds_zm(2)========== m(2)
  //
  //    --xm(1)---invrs_rho_ds_zt(1)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(1)
  //
  //    ========[x'w'|_sfc]=====================rho_ds_zm(1)========== m(1) sfc
  //
  //    where "sfc" is the level of the model surface or lower boundary.
  //
  //    The vertically discretized form of the turbulent advection term
  //    (treated as an eddy diffusion term), with the explicit surface flux,
  //    x'w'|_sfc, in place, is written out as:
  //
  //    - invrs_rho_ds_zt(1)
  //      * invrs_dzt(1) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];
  //
  //    and after substituting x'w'(2) = -K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) ),
  //    this expression can be re-written as:
  //
  //    + invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [   rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
  //            + rho_ds_zm(1) * x'w'|_sfc ];
  //
  //    which can be re-written again as:
  //
  //    + invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
  //    + invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(1) * x'w'|_sfc.
  //
  //    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
  //    is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
  //    eddy-diffusion term.  The discretized implicit portion of the term is
  //    written out as:
  //
  //    + (1/2) * invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [ rho_ds_zm(2) * K_zm(2)
  //            * invrs_dzm(2) * ( xm(2,<t+1>) - xm(1,<t+1>) ) ].
  //
  //    Note:  When the implicit term is brought over to the left-hand side,
  //           the sign is reversed and the leading "+" in front of the term
  //           is changed to a "-".
  //
  //    The discretized explicit portion of the term is written out as:
  //
  //    + (1/2) * invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [ rho_ds_zm(2) * K_zm(2)
  //            * invrs_dzm(2) * ( xm(2,<t>) - xm(1,<t>) ) ]
  //    + invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(1) * x'w'|_sfc.
  //
  //    Note:  The x'w'|_sfc portion of the term written above has been pulled
  //           away from the rest of the explicit form written above because
  //           the (1/2) factor due to Crank-Nicholson time_stepping does not
  //           apply to it, as there isn't an implicit portion for x'w'|_sfc.
  //
  //    Timestep index (t) stands for the index of the current timestep, while
  //    timestep index (t+1) stands for the index of the next timestep, which
  //    is being advanced to in solving the d(xm)/dt equation.
  //
  // 2) x'w' surface flux; implicit form for momentum fluxes u'w' and v'w'
  //
  //    The x'w' surface flux is applied to the d(xm)/dt equation through the
  //    turbulent advection term, which is:
  //
  //    - (1/rho_ds) * d( rho_ds * x'w' )/dz.
  //
  //    At most vertical levels, a substitution can be made for x'w', such
  //    that:
  //
  //    x'w' = - K_zm * d(xm)/dz.
  //
  //    However, the same substitution cannot be made at the surface (momentum
  //    level 1), as x'w'|_sfc is a surface momentum flux that is found by the
  //    following equation:
  //
  //    x'w'|_sfc = - [ u_star^2 / sqrt( um^2 + vm^2 ) ] * xm;
  //
  //    where x'w'|_sfc and xm are either u'w'|_sfc and um, respectively, or
  //    v'w'|_sfc and vm, respectively (um and vm are located at the first
  //    thermodynamic level above the surface, which is thermodynamic level 1),
  //    sqrt( um^2 + vm^2 ) is the wind speed (also at thermodynamic level 1),
  //    and u_star is defined as:
  //
  //    u_star = ( u'w'|_sfc^2 + v'w'|_sfc^2 )^(1/4);
  //
  //    and thus u_star^2 is defined as:
  //
  //    u_star^2 = sqrt( u'w'|_sfc^2 + v'w'|_sfc^2 ).
  //
  //    The value of u_star is either set to a constant value or computed
  //    (through function diag_ustar) based on the surface wind speed, the
  //    height above surface of the surface wind speed (as compared to the
  //    roughness height), and the buoyancy flux at the surface.  Either way,
  //    u_star is computed elsewhere in the model, and the values of u'w'|_sfc
  //    and v'w'|_sfc are based on it and computed along with it.  The values
  //    of u'w'|_sfc and v'w'|_sfc are then passed into advance_clubb_core,
  //    and are eventually passed into advance_windm_edsclrm.  In subroutine
  //    advance_windm_edsclrm, the value of u_star_sqd is then recomputed
  //    based on u'w'|_sfc and v'w'|_sfc.  The value of sqrt( u_star_sqd ) is
  //    consistent with the value of the original computation of u_star.
  //
  //    The equation listed above is substituted for x'w'|_sfc.  The lower
  //    boundary condition, which in this case needs to be applied to the
  //    d(xm)/dt equation at level 1, is discretized as follows:
  //
  //    --xm(2)------------------------------------------------------- t(2)
  //
  //    ===[x'w'(2) = -K_zm(2)*d(xm)/dz]=================rho_ds_zm(2)= m(2)
  //
  //    --xm(1)---invrs_rho_ds_zt(1)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(1)
  //
  //    ===[x'w'|_sfc = -[u_star^2/sqrt(um^2+vm^2)]*xm]==rho_ds_zm(1)= m(1) sfc
  //
  //    where "sfc" is the level of the model surface or lower boundary.
  //
  //    The vertically discretized form of the turbulent advection term
  //    (treated as an eddy diffusion term), with the implicit surface momentum
  //    flux in place, is written out as:
  //
  //    - invrs_rho_ds_zt(1)
  //      * invrs_dzt(1) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];
  //
  //    which can be re-written as:
  //
  //    - invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [   rho_ds_zm(2)
  //              * { - K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) ) }
  //            - rho_ds_zm(1)
  //              * { - [ u_star^2 / sqrt( um(1)^2 + vm(1)^2 ) ] * xm(1) } ];
  //
  //    which can be re-written as:
  //
  //    + invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
  //    - invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(1) * [ u_star^2 / sqrt( um(1)^2 + vm(1)^2 ) ] * xm(1).
  //
  //    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
  //    is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
  //    eddy-diffusion term.  The discretized implicit portion of the term is
  //    written out as:
  //
  //    + (1/2) * invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [ rho_ds_zm(2) * K_zm(2)
  //            * invrs_dzm(2) * ( xm(2,<t+1>) - xm(1,<t+1>) ) ]
  //    - invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * rho_ds_zm(1)
  //        * [u_star^2/sqrt( um(1,<t>)^2 + vm(1,<t>)^2 )] * xm(1,<t+1>).
  //
  //    Note:  When the implicit term is brought over to the left-hand side,
  //           the signs are reversed and the leading "+" in front of the first
  //           part of the term is changed to a "-", while the leading "-" in
  //           front of the second part of the term is changed to a "+".
  //
  //    Note:  The x'w'|_sfc portion of the term written above has been pulled
  //           away from the rest of the implicit form written above because
  //           the (1/2) factor due to Crank-Nicholson time_stepping does not
  //           apply to it.  The x'w'|_sfc portion of the term is treated
  //           completely implicitly in order to enhance numerical stability.
  //
  //    The discretized explicit portion of the term is written out as:
  //
  //    + (1/2) * invrs_rho_ds_zt(1)
  //      * invrs_dzt(1)
  //        * [ rho_ds_zm(2) * K_zm(2)
  //            * invrs_dzm(2) * ( xm(2,<t>) - xm(1,<t>) ) ].
  //
  //    Timestep index (t) stands for the index of the current timestep, while
  //    timestep index (t+1) stands for the index of the next timestep, which
  //    is being advanced to in solving the d(xm)/dt equation.
  //
  //
  // The lower boundary condition for the implicit and explicit portions of the
  // turbulent advection term, without the x'w'|_sfc portion of the term, can
  // easily be invoked by using the zero-flux boundary conditions found in the
  // generalized diffusion function (function diffusion_zt_lhs), which is used
  // for many other equations in this model.  Either the generalized explicit
  // surface flux needs to be added onto the explicit term after the diffusion
  // function has been called from subroutine windm_edsclrm_rhs, or the
  // implicit momentum surface flux needs to be added onto the implicit term
  // after the diffusion function has been called from subroutine
  // windm_edsclrm_lhs.
  //
  //
  // Conservation Properties:
  //
  // When a fixed-flux lower boundary condition is used (combined with a
  // zero-flux upper boundary condition), this technique of discretizing the
  // turbulent advection term (treated as an eddy-diffusion term) leads to
  // conservative differencing.  When the implicit momentum surface flux is
  // either zero or not used, the column totals for each column in the
  // left-hand side matrix (for the turbulent advection term) should be equal
  // to 0.  Otherwise, the column total for the second column will be equal to
  // rho_ds_zm(1) * x'w'|_sfc<t+1>.   When the generalized explicit surface
  // flux is either zero or not used, the column total for the right-hand side
  // vector (for the turbulent advection term) should be equal to 0.
  // Otherwise, the column total for the right-hand side vector (for the
  // turbulent advection term) will be equal to rho_ds_zm(1) * x'w'|_sfc<t>.
  // This ensures that the total amount of quantity xm over the entire vertical
  // domain is only changed by the surface flux (neglecting any forcing terms).
  // The total amount of change is equal to rho_ds_zm(1) * x'w'|_sfc.
  //
  // To see that this conservation law is satisfied by the left-hand side
  // matrix, compute the turbulent advection (treated as eddy diffusion) of xm,
  // neglecting any implicit momentum surface flux, multiply by rho_ds_zt, and
  // integrate vertically.  In discretized matrix notation (where "i" stands
  // for the matrix column and "j" stands for the matrix row):
  //
  //  0 = Sum_j Sum_i
  //       (rho_ds_zt)_i ( dzt )_i
  //       ( 0.5 * (1/rho_ds_zt) * (1/dzt) * (rho_ds_zm*K_zm*(1/dzm)) )_ij (xm<t+1>)_j.
  //
  // The left-hand side matrix,
  // ( 0.5 * (1/rho_ds_zt) * (1/dzt) * (rho_ds_zm*K_zm*(1/dzm)) )_ij, is
  // partially written below.  The sum over i in the above equation removes
  // (1/rho_ds_zt) and (1/dzt) everywhere from the matrix below.  The sum
  // over j leaves the column totals that are desired, which are 0.
  //
  // Left-hand side matrix contributions from the turbulent advection term
  // (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
  // first four vertical levels:
  //
  //     -------------------------------------------------------------------------------->
  // k=1 |  +0.5*                    -0.5*                                0
  //    |    (1/rho_ds_zt(k))*        (1/rho_ds_zt(k))*
  //    |    (1/dzt(k))*              (1/dzt(k))*
  //    |    rho_ds_zm(k+1)*          rho_ds_zm(k+1)*
  //    |    K_zm(k+1)*(1/dzm(k+1))   K_zm(k+1)*(1/dzm(k+1))
  //    |
  // k=2 |  -0.5*                    +0.5*                      -0.5*
  //    |    (1/rho_ds_zt(k))*        (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
  //    |    (1/dzt(k))*              (1/dzt(k))*                (1/dzt(k))*
  //    |    rho_ds_zm(k)*            [ rho_ds_zm(k+1)*          rho_ds_zm(k+1)*
  //    |    K_zm(k)*(1/dzm(k))         K_zm(k+1)*(1/dzm(k+1))   K_zm(k+1)*(1/dzm(k+1))
  //    |                              +rho_ds_zm(k)*
  //    |                               K_zm(k)*(1/dzm(k)) ]
  //    |
  // k=3 |            0              -0.5*                      +0.5*
  //    |                             (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
  //    |                             (1/dzt(k))*                (1/dzt(k))*
  //    |                             rho_ds_zm(k)*              [ rho_ds_zm(k+1)*
  //    |                             K_zm(k)*(1/dzm(k))           K_zm(k+1)*(1/dzm(k+1))
  //    |                                                         +rho_ds_zm(k)*
  //    |                                                          K_zm(k)*(1/dzm(k)) ]
  //    |
  // k=4 |            0                        0                -0.5*
  //    |                                                        (1/rho_ds_zt(k))*
  //    |                                                        (1/dzt(k))*
  //    |                                                        rho_ds_zm(k)*
  //    |                                                        K_zm(k)*(1/dzm(k))
  //   \ /
  //
  // Note:  The superdiagonal term from level 3 and both the main diagonal and
  //        superdiagonal terms from level 4 are not shown on this diagram.
  //
  // Note:  If an implicit momentum surface flux is used, an additional term,
  //        + (1/rho_ds_zt(1)) * (1/dzt(1)) * rho_ds_zm(1)
  //          * [ u_star^2 / sqrt( um(1,<t>)^2 + vm(1,<t>)^2 ) ], is added to
  //        row 1 (k=1), column 1.
  //
  // To see that the above conservation law is satisfied by the right-hand side
  // vector, compute the turbulent advection (treated as eddy diffusion) of xm,
  // neglecting any generalized explicit surface flux, multiply by rho_ds_zt,
  // and integrate vertically.  In discretized matrix notation (where "i"
  // stands for the matrix column and "j" stands for the matrix row):
  //
  //  0 = Sum_j Sum_i (rho_ds_zt)_i ( dzt )_i ( rhs_vector )_j.
  //
  // The right-hand side vector, ( rhs_vector )_j, is partially written below.
  // The sum over i in the above equation removes (1/rho_ds_zt) and (1/dzt)
  // everywhere from the vector below.  The sum over j leaves the column total
  // that is desired, which is 0.
  //
  // Right-hand side vector contributions from the turbulent advection term
  // (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
  // first four vertical levels:
  //
  //     --------------------------------------------------
  // k=1 | +0.5*(1/rho_ds_zt(k))*                           |
  //    |      (1/dzt(k))*                                 |
  //    |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
  //    |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>)) ]   |
  //    |                                                  |
  // k=2 | +0.5*(1/rho_ds_zt(k))*                           |
  //    |      (1/dzt(k))*                                 |
  //    |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
  //    |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
  //    |        -rho_ds_zm(k)*K_zm(k)*                    |
  //    |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
  //    |                                                  |
  // k=3 | +0.5*(1/rho_ds_zt(k))*                           |
  //    |      (1/dzt(k))*                                 |
  //    |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
  //    |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
  //    |        -rho_ds_zm(k)*K_zm(k)*                    |
  //    |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
  //    |                                                  |
  // k=4 | +0.5*(1/rho_ds_zt(k))*                           |
  //    |      (1/dzt(k))*                                 |
  //    |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
  //    |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
  //    |        -rho_ds_zm(k)*K_zm(k)*                    |
  //    |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
  //   \ /                                                \ /
  //
  // Note:  If a generalized explicit surface flux is used, an additional term,
  //        + (1/rho_ds_zt(1)) * (1/dzt(1)) * rho_ds_zm(1) * x'w'|_sfc, is
  //        added to row 1 (k=1).
  //
  // Note:  Only the contributions by the turbulent advection term are shown
  //        for both the left-hand side matrix and the right-hand side vector.
  //        There are more terms in the equation, and thus more factors to be
  //        added to both the left-hand side matrix (such as time tendency and
  //        mean advection) and the right-hand side vector (such as xm
  //        forcings).  The left-hand side matrix is set-up so that a singular
  //        matrix is not encountered.
  // References:
  // Eqn. 8 & 9 on p. 3545 of
  // ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  //   Method and Model Description'' Golaz, et al. (2002)
  // JAS, Vol. 59, pp. 3540--3551.
  //-----------------------------------------------------------------------
  // ------------------------ Input Variables ------------------------
  // ------------------------ Inout variables ------------------------
  // ------------------------ Output variables ------------------------
  // ------------------------ Local variables ------------------------
  // Est. of the condition number of the variance LHS matrix
  auto workspace_scope = workspace.scope();
  Kokkos::View<real_t*, Kokkos::LayoutLeft> rcond =
      workspace_scope.view<Kokkos::View<real_t*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "windm_edsclrm_solve::rcond",
          ngrdcol); // Estimate of the reciprocal of the condition number on the LHS matrix

  // ------------------------ Begin Code ------------------------

  // Matrix solves are bit-different between ascending and descending.
  // This ensures matrices are solved in the same (descending) order,
  // which is useful for ensuring BFBness between grid modes
  // Native lhs is (band,column,row); transpose only at the existing
  // (column,band,row) matrix API boundary. Preserve Native descending flips.
  const bool reverse = l_force_descending_solves && gr.grid_dir_indx > 0;
  Real3D lhs_solver  = workspace_scope.view<Real3D>(
      Kokkos::WithoutInitializing, "windm_edsclrm_solve::lhs_solver", ngrdcol, 3, nzt);
  Real3D rhs_solver = workspace_scope.view<Real3D>(
      Kokkos::WithoutInitializing, "windm_edsclrm_solve::rhs_solver", ngrdcol, nzt, nrhs);
  Real3D solution_solver = workspace_scope.view<Real3D>(
      Kokkos::WithoutInitializing, "windm_edsclrm_solve::solution_solver", ngrdcol, nzt, nrhs);
  Kokkos::parallel_for(
      "windm_edsclrm_solve::matrix_order", (ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
        const int i = (cell % ngrdcol), k = (cell / ngrdcol),
                  source_k = reverse ? (-k + nzt - 1) : k;
        for (int b = 0; b < 3; ++b)
          lhs_solver(i, b, k) = lhs(reverse ? (-b + 2) : b, i, source_k);
        for (int n = 0; n < nrhs; ++n)
          rhs_solver(i, k, n) = rhs(i, source_k, n);
      });

  if (l_need_rcond) {

    tridiag_solve(workspace, "windm_edsclrm", tridiag_solve_method, ngrdcol, nzt, nrhs,
                  l_implemented,                    // In
                  lhs_solver, rhs_solver, err_info, // InOut
                  solution_solver,                  // Out
                  &rcond);                          // Optional Out
    // Copy InOut allocations and the unflipped answer before Native diagnostics.
    Kokkos::parallel_for(
        "windm_edsclrm_solve::native_storage", (ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          for (int b = 0; b < 3; ++b)
            lhs(b, i, k) = lhs_solver(i, b, k);
          for (int n = 0; n < nrhs; ++n) {
            rhs(i, k, n)      = rhs_solver(i, k, n);
            solution(i, k, n) = solution_solver(i, k, n);
          }
        });

    {
      Real1D stats_values_1 = workspace_scope.view<Real1D>(
          Kokkos::WithoutInitializing, "windm_edsclrm_solve::stats_values_1", rcond.extent(0));
      Kokkos::parallel_for(
          "windm_edsclrm_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
            const int i       = cell;
            stats_values_1(i) = real_t(1.0) / rcond(i);
          });
      stats_update("windm_matrix_condt_num", stats_values_1, // In
                   stats);                                   // InOut
    }
  } else {

    tridiag_solve(workspace, "windm_edsclrm", tridiag_solve_method, ngrdcol, nzt, nrhs,
                  l_implemented,                    // In
                  lhs_solver, rhs_solver, err_info, // InOut
                  solution_solver,                  // Out
                  nullptr);                         // Optional Out
    // Copy InOut allocations and the unflipped answer before Native diagnostics.
    Kokkos::parallel_for(
        "windm_edsclrm_solve::native_storage", (ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          for (int b = 0; b < 3; ++b)
            lhs(b, i, k) = lhs_solver(i, b, k);
          for (int n = 0; n < nrhs; ++n) {
            rhs(i, k, n)      = rhs_solver(i, k, n);
            solution(i, k, n) = solution_solver(i, k, n);
          }
        });
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling tridiag_solve in windm_edsclrm_solve" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  // Flip the back to the ascending direction if we forced the solve
  // to be in descending mode
  // Restore direction only after Native debug-fatal return handling.
  if (l_force_descending_solves && gr.grid_dir_indx > 0) {
    Kokkos::parallel_for(
        "windm_edsclrm_solve::restore_direction", (ngrdcol * nzt), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          for (int n = 0; n < nrhs; ++n)
            solution(i, k, n) = solution_solver(i, (-k + nzt - 1), n);
        });
  }

  return;
}

//=============================================================================
// solve_type: Desc. of what is being solved for
// xm: Computed value um or vm at <t+1>    [m/s]
// invrs_dzt: The inverse spacing between momentum grid levels; centered over thermodynamic grid
// levels.
// lhs_diff: LHS diffustion terms
// lhs_ma_zt: LHS mean advection terms
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels  [m^3/kg]
// u_star_sqd: Surface friction velocity, u_star, squared      [m^2/s^2]
// rho_ds_zm: Dry, static density on momentum levels      [kg/m^3]
// wind_speed: wind speed; sqrt(u^2 + v^2)              [m/s]
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7>
void windm_edsclrm_implicit_stats(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                  int solve_type,
                                  const grid& gr,                            // In
                                  Array0 xm, Array1 invrs_dzt,               // In
                                  Array2 lhs_diff, Array3 lhs_ma_zt,         // In
                                  Array4 invrs_rho_ds_zt, Array5 u_star_sqd, // In
                                  Array6 rho_ds_zm, Array7 wind_speed,       // In
                                  bool l_imp_sfc_momentum_flux,              // In
                                  stats_type& stats)                         // InOut
{
  // Description:
  // Compute implicit contributions to um and vm
  // References:
  // None
  //-----------------------------------------------------------------------
  //---------------------- Input Variables ----------------------
  //---------------------- Local variables ----------------------
  // Budget names
  // km1: Array indices
  auto workspace_scope = workspace.scope();
  std::string name_ma;
  std::string name_ta;
  const auto gr_k_lb_zt       = gr.k_lb_zt;
  const auto gr_grid_dir      = gr.grid_dir;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_ub_zt       = gr.k_ub_zt;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> imp_sfc_flux =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "windm_edsclrm_implicit_stats::imp_sfc_flux", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "windm_edsclrm_implicit_stats::stats_tmp", ngrdcol, nzt);

  //---------------------- Begin Code ----------------------

  if (!stats.l_sample) {
    return;
  }

  switch (solve_type) {
  case windm_edsclrm_um: {
    name_ma = "um_ma";
    name_ta = "um_ta";

    break;
  }
  case windm_edsclrm_vm: {
    name_ma = "vm_ma";
    name_ta = "vm_ta";

    break;
  }
  default: {
    name_ma = "";
    name_ta = "";
    return;
  }
  }

  Kokkos::parallel_for(
      "windm_edsclrm_implicit_stats::array1", (imp_sfc_flux.extent(0) * imp_sfc_flux.extent(1)),
      KOKKOS_LAMBDA(int cell) {
        const int slice_0              = ((cell % imp_sfc_flux.extent(0)));
        const int slice_1              = ((cell / imp_sfc_flux.extent(0)));
        imp_sfc_flux(slice_0, slice_1) = zero;
      });

  if (l_imp_sfc_momentum_flux) {

    // Statistics:  implicit contributions for um or vm.
    // xm term ta is modified at level 2 to include the effects of the
    // surface flux.  In this case, this affects the implicit portion of
    // the term, which handles the main diagonal for the turbulent advection term.
    if (var_on_stats_list(stats, name_ta)) {
      Kokkos::parallel_for(
          "windm_edsclrm_implicit_stats::loop2", std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                 = cell;
            imp_sfc_flux(i, gr_k_lb_zt) = -gr_grid_dir * invrs_rho_ds_zt(i, gr_k_lb_zt) *
                                          invrs_dzt(i, gr_k_lb_zt) * rho_ds_zm(i, gr_k_lb_zm) *
                                          (u_star_sqd(i) / wind_speed(i, gr_k_lb_zt));
          });
    }

    // l_imp_sfc_momentum_flux
  }

  if (stats.l_sample) {
    Kokkos::parallel_for(
        "windm_edsclrm_implicit_stats::array3", (stats_tmp.extent(0) * stats_tmp.extent(1)),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0           = ((cell % stats_tmp.extent(0)));
          const int slice_1           = ((cell / stats_tmp.extent(0)));
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "windm_edsclrm_implicit_stats::loop4", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          int k;
          int kp1;
          int km1;
          // Lower boundary conditions
          // xm mean advection
          // xm term ma is completely implicit; record with stats_update.
          k               = gr_k_lb_zt;
          kp1             = (gr_grid_dir_indx + k);
          stats_tmp(i, k) = -lhs_ma_zt(1, i, k) * xm(i, k) -
                            lhs_ma_zt((-gr_grid_dir_indx + 1), i, k) * xm(i, kp1);

          // Interior domain
          // xm mean advection
          // xm term ma is completely implicit; record with stats_update.
          for (k = 1; k < (nzt - 1); k += 1) {
            km1             = (-gr_grid_dir_indx + k);
            kp1             = (gr_grid_dir_indx + k);
            stats_tmp(i, k) = -lhs_ma_zt((gr_grid_dir_indx + 1), i, k) * xm(i, km1) -
                              lhs_ma_zt(1, i, k) * xm(i, k) -
                              lhs_ma_zt((-gr_grid_dir_indx + 1), i, k) * xm(i, kp1);
          }

          // Upper boundary conditions
          // xm mean advection
          // xm term ma is completely implicit; record with stats_update.
          k   = gr_k_ub_zt;
          km1 = (-gr_grid_dir_indx + k);
          stats_tmp(i, k) =
              -lhs_ma_zt((gr_grid_dir_indx + 1), i, k) * xm(i, km1) - lhs_ma_zt(1, i, k) * xm(i, k);
        });
    stats_update(name_ma, stats_tmp, // In
                 stats);             // InOut

    Kokkos::parallel_for(
        "windm_edsclrm_implicit_stats::array5", (stats_tmp.extent(0) * stats_tmp.extent(1)),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0           = ((cell % stats_tmp.extent(0)));
          const int slice_1           = ((cell / stats_tmp.extent(0)));
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "windm_edsclrm_implicit_stats::loop6", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          int k;
          int kp1;
          int km1;
          // Lower boundary conditions
          // xm turbulent transport (implicit component)
          // xm term ta has both implicit and explicit components;
          // finalize the implicit part with stats_finalize_budget.
          k               = gr_k_lb_zt;
          kp1             = (gr_grid_dir_indx + k);
          stats_tmp(i, k) = (-real_t(0.5) * lhs_diff(1, i, k) + imp_sfc_flux(i, k)) * xm(i, k) -
                            real_t(0.5) * lhs_diff((-gr_grid_dir_indx + 1), i, k) * xm(i, kp1);

          // Interior domain
          // xm turbulent transport (implicit component)
          // xm term ta has both implicit and explicit components;
          // finalize the implicit part with stats_finalize_budget.
          for (k = 1; k < (nzt - 1); k += 1) {
            km1             = (-gr_grid_dir_indx + k);
            kp1             = (gr_grid_dir_indx + k);
            stats_tmp(i, k) = -real_t(0.5) * lhs_diff((gr_grid_dir_indx + 1), i, k) * xm(i, km1) +
                              (-real_t(0.5) * lhs_diff(1, i, k) + imp_sfc_flux(i, k)) * xm(i, k) -
                              real_t(0.5) * lhs_diff((-gr_grid_dir_indx + 1), i, k) * xm(i, kp1);
          }

          // Upper boundary conditions
          // xm turbulent transport (implicit component)
          // xm term ta has both implicit and explicit components;
          // finalize the implicit part with stats_finalize_budget.
          k               = gr_k_ub_zt;
          km1             = (-gr_grid_dir_indx + k);
          stats_tmp(i, k) = -real_t(0.5) * lhs_diff((gr_grid_dir_indx + 1), i, k) * xm(i, km1) +
                            (-real_t(0.5) * lhs_diff(1, i, k) + imp_sfc_flux(i, k)) * xm(i, k);
        });
    // Finalize implicit contributions for xm
    stats_finalize_budget(name_ta, stats_tmp, // In
                          stats);             // InOut
  }

  return;
}

//=============================================================================
// solve_type: Description of what is being solved for
// fcor: Coriolis parameter     [s^-1]
// perp_wind_m: Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
// perp_wind_g: Perpendicular component of the geostropic wind (e.g. vg)         [m/s]
// xm_forcing: Prescribed wind forcing                                          [m/s/s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// xm_tndcy: xm tendency            [m/s^2]
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void compute_uv_tndcy(ScratchWorkspace& workspace,                         // InOut
                      int nzt, int ngrdcol, int solve_type,                // In
                      Array0 fcor, Array1 perp_wind_m, Array2 perp_wind_g, // In
                      Array3 xm_forcing, bool l_implemented,               // In
                      stats_type& stats,                                   // InOut
                      Array4 xm_tndcy)                                     // Out
{
  // Description:
  // Computes the explicit tendency for the um and vm wind components.
  //
  // The only explicit tendency that is involved in the d(um)/dt or d(vm)/dt
  // equations is the Coriolis tendency.
  //
  // The d(um)/dt equation contains the term:
  //
  // - f * ( v_g - vm );
  //
  // where f is the Coriolis parameter and v_g is the v component of the
  // geostrophic wind.
  //
  // Likewise, the d(vm)/dt equation contains the term:
  //
  // + f * ( u_g - um );
  //
  // where u_g is the u component of the geostrophic wind.
  //
  // This term is treated completely explicitly.  The values of um, vm, u_g,
  // and v_g are all found on the thermodynamic levels.
  //
  // Wind forcing from the GCSS cases is also added here.
  //
  // References:
  //-----------------------------------------------------------------------
  // -------------------------- Input Variables --------------------------
  // -------------------------- Output Variables --------------------------
  // -------------------------- Local Variables --------------------------
  auto workspace_scope = workspace.scope();
  std::string name_gf;
  std::string name_cf;
  std::string name_f;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_gf =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "compute_uv_tndcy::xm_gf", ngrdcol, nzt);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_cf =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "compute_uv_tndcy::xm_cf", ngrdcol, nzt);

  // -------------------------- Begin Code --------------------------

  if (!l_implemented) {
    // Only compute the Coriolis term if the model is running on it's own,
    // and is not part of a larger, host model.

    switch (solve_type) {
    case windm_edsclrm_um: {

      name_gf = "um_gf";
      name_cf = "um_cf";
      name_f  = "um_f";

      Kokkos::parallel_for(
          "compute_uv_tndcy::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            xm_gf(i, k) = -fcor(i) * perp_wind_g(i, k);
          });

      Kokkos::parallel_for(
          "compute_uv_tndcy::loop2", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            xm_cf(i, k) = fcor(i) * perp_wind_m(i, k);
          });

      break;
    }
    case windm_edsclrm_vm: {

      name_gf = "vm_gf";
      name_cf = "vm_cf";
      name_f  = "vm_f";

      Kokkos::parallel_for(
          "compute_uv_tndcy::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            xm_gf(i, k) = fcor(i) * perp_wind_g(i, k);
          });

      Kokkos::parallel_for(
          "compute_uv_tndcy::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            xm_cf(i, k) = -fcor(i) * perp_wind_m(i, k);
          });

      break;
    }
    default: {

      name_gf = "";
      name_cf = "";
      name_f  = "";

      Kokkos::parallel_for(
          "compute_uv_tndcy::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = ((cell % ngrdcol));
            const int k = ((cell / ngrdcol));
            xm_gf(i, k) = real_t(0.);
            xm_cf(i, k) = real_t(0.);
          });

      break;
    }
    }

    Kokkos::parallel_for(
        "compute_uv_tndcy::loop6", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i    = ((cell % ngrdcol));
          const int k    = ((cell / ngrdcol));
          xm_tndcy(i, k) = xm_gf(i, k) + xm_cf(i, k) + xm_forcing(i, k);
        });

    if (stats.l_sample) {
      // Statistics:  explicit contributions for um or vm.
      // xm term gf is completely explicit; record with stats_update.
      stats_update(name_gf, xm_gf, // In
                   stats);         // InOut
      // xm term cf is completely explicit; record with stats_update.
      stats_update(name_cf, xm_cf, // In
                   stats);         // InOut
      // xm term F
      stats_update(name_f, xm_forcing, // In
                   stats);             // InOut
    }

  } else {
    // implemented in a host model.

    Kokkos::parallel_for(
        "compute_uv_tndcy::loop7", std::max(0, nzt) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i    = ((cell % ngrdcol));
          const int k    = ((cell / ngrdcol));
          xm_tndcy(i, k) = real_t(0.0);
        });
  }

  return;
}

//=============================================================================
// dt: Model timestep                             [s]
// lhs_ma_zt: LHS mean advection terms
// lhs_diff: LHS diffustion terms
// wind_speed: wind speed; sqrt( u^2 + v^2 )              [m/s]
// u_star_sqd: Surface friction velocity, u_*, squared  [m^2/s^2]
// rho_ds_zm: Dry, static density on momentum levels     [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels [m^3/kg]
// l_implemented: True if CLUBB is being implemented and run in a host model
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
// lhs: Implicit contributions to xm (tridiagonal matrix)
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void windm_edsclrm_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, // In
                       Array0 lhs_ma_zt, Array1 lhs_diff,                        // In
                       Array2 wind_speed, Array3 u_star_sqd,                     // In
                       Array4 rho_ds_zm, Array5 invrs_rho_ds_zt,                 // In
                       bool l_implemented, bool l_imp_sfc_momentum_flux,         // In
                       Array6 lhs)                                               // Out
{
  // Description:
  // Calculate the implicit portion of the horizontal wind or eddy-scalar
  // time-tendency equation.  See the description in subroutine
  // windm_edsclrm_solve for more details.
  //
  //   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  //   Simple changes to this procedure may adversely affect computational speed
  //       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  //----------------------------------------------------------------------------------
  // ----------------------- Input Variables -----------------------
  // ----------------------- Output Variable -----------------------
  // ----------------------- Local Variables -----------------------
  real_t invrs_dt; // Inverse of dt, 1/dt, used for computational efficiency
  // i: Loop variable
  const auto gr_k_lb_zt       = gr.k_lb_zt;
  const auto gr_k_ub_zt       = gr.k_ub_zt;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_grid_dir      = gr.grid_dir;
  const auto gr_invrs_dzt     = gr.invrs_dzt;
  const auto gr_k_lb_zm       = gr.k_lb_zm;

  // ----------------------- Begin Code -----------------------

  // Calculate coefs of eddy diffusivity and inverse of dt
  invrs_dt = real_t(1.0) / dt;

  // Add terms to lhs
  Kokkos::parallel_for(
      "windm_edsclrm_lhs::loop1", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i  = ((cell % ngrdcol));
        const int k  = ((cell / ngrdcol));
        lhs(0, i, k) = real_t(0.5) * lhs_diff(0, i, k);

        lhs(1, i, k) = real_t(0.5) * lhs_diff(1, i, k);

        // LHS time tendency.
        lhs(1, i, k) = lhs(1, i, k) + invrs_dt;

        lhs(2, i, k) = real_t(0.5) * lhs_diff(2, i, k);
      });

  // LHS mean advection term.
  if (!l_implemented) {
    Kokkos::parallel_for(
        "windm_edsclrm_lhs::loop2",
        std::max(0, ((((-gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1))) *
            std::max(0, ngrdcol) * std::max(0, 3),
        KOKKOS_LAMBDA(const int cell) {
          const int j  = ((cell % 3));
          const int i  = (((cell / 3) % ngrdcol));
          const int k  = (gr_grid_dir_indx * (cell / (3 * ngrdcol)) + gr_k_lb_zt);
          lhs(j, i, k) = lhs(j, i, k) + lhs_ma_zt(j, i, k);
        });
  }

  if (l_imp_sfc_momentum_flux) {

    // LHS momentum surface flux.
    Kokkos::parallel_for(
        "windm_edsclrm_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          lhs(1, i, gr_k_lb_zt) =
              lhs(1, i, gr_k_lb_zt) + gr_grid_dir * invrs_rho_ds_zt(i, gr_k_lb_zt) *
                                          gr_invrs_dzt(i, gr_k_lb_zt) * rho_ds_zm(i, gr_k_lb_zm) *
                                          (u_star_sqd(i) / wind_speed(i, gr_k_lb_zt));
        });

    // l_imp_sfc_momentum_flux
  }

  return;
}

//=============================================================================
// solve_type: Description of what is being solved for
// dt: Model timestep                             [s]
// lhs_diff: LHS diffusion terms
// xm: Eddy-scalar variable, xm (thermo. levels)  [units vary]
// xm_tndcy: The explicit time-tendency acting on xm    [units vary]
// rho_ds_zm: Dry, static density on momentum levels     [kg/m^3]
// invrs_rho_ds_zt: Inv. dry, static density at thermo. levels [m^3/kg]
// l_imp_sfc_momentum_flux: Flag for implicit momentum surface fluxes.
// xpwp_sfc: x'w' at the surface                              [units vary]
// rhs: Right-hand side (explicit) contributions.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void windm_edsclrm_rhs(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                       int solve_type,
                       real_t dt,                                     // In
                       Array0 lhs_diff, Array1 xm, Array2 xm_tndcy,   // In
                       Array3 rho_ds_zm, Array4 invrs_rho_ds_zt,      // In
                       bool l_imp_sfc_momentum_flux, Array5 xpwp_sfc, // In
                       stats_type& stats,                             // InOut
                       Array6 rhs)                                    // Out
{
  // Description:
  //   Calculate the explicit portion of the horizontal wind or eddy-scalar
  //   time-tendency equation.  See the description in subroutine
  //   windm_edsclrm_solve for more details.
  //
  // References:
  //   None
  //
  //   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  //   Simple changes to this procedure may adversely affect computational speed
  //       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  //------------------- Inout Variable -------------------
  //------------------- Output Variable -------------------
  //------------------- Local Variables -------------------
  // k: Loop variable
  auto workspace_scope = workspace.scope();
  real_t invrs_dt;
  std::string name_ta;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zt       = gr.k_lb_zt;
  const auto gr_grid_dir      = gr.grid_dir;
  const auto gr_invrs_dzt     = gr.invrs_dzt;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp =
      workspace_scope.view<Kokkos::View<real_t**, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "windm_edsclrm_rhs::stats_tmp", ngrdcol, nzt);

  //------------------- Begin Code -------------------

  // Precalculate 1.0/dt to avoid redoing the divide
  invrs_dt = real_t(1.0) / dt;

  switch (solve_type) {
  case windm_edsclrm_um: {
    name_ta = "um_ta";
    break;
  }
  case windm_edsclrm_vm: {
    name_ta = "vm_ta";
    break;
  }
  default: {
    // Eddy scalars
    name_ta = "";
    break;
  }
  }

  // Lower (ascending grid) or upper (descending grid) boundary calculation
  Kokkos::parallel_for(
      "windm_edsclrm_rhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        // RHS forcings
        // RHS time tendency
        rhs(i, 0) = real_t(0.5) * (-lhs_diff(1, i, 0) * xm(i, 0) - lhs_diff(0, i, 0) * xm(i, 1)) +
                    xm_tndcy(i, 0) + invrs_dt * xm(i, 0);
      });

  // Non-boundary rhs calculation, this is a highly vectorized loop
  // Note: Some gymnastics are done here to produce a bit-for-bit match
  //       in the ascending vs. descending grid test at -O0 optimization.
  //       This code forces the same calculations to be done in the same
  //       order, regardless of grid direction.
  // For ascending:
  // ( - lhs_diff(3,i,k) * xm(i,k-1)
  //   - lhs_diff(2,i,k) * xm(i,k)
  //   - lhs_diff(1,i,k) * xm(i,k+1) )
  // For descending:
  // ( - lhs_diff(1,i,k) * xm(i,k+1)
  //   - lhs_diff(2,i,k) * xm(i,k)
  //   - lhs_diff(3,i,k) * xm(i,k-1) )
  Kokkos::parallel_for(
      "windm_edsclrm_rhs::loop2", std::max(0, ((nzt - 2))) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        const int i = ((cell % ngrdcol));
        const int k = ((cell / ngrdcol) + 1);
        // RHS forcings
        // RHS time tendency
        rhs(i, k) = real_t(0.5) *
                        (-lhs_diff((gr_grid_dir_indx + 1), i, k) * xm(i, (-gr_grid_dir_indx + k)) -
                         lhs_diff(1, i, k) * xm(i, k) -
                         lhs_diff((-gr_grid_dir_indx + 1), i, k) * xm(i, (gr_grid_dir_indx + k))) +
                    xm_tndcy(i, k) + invrs_dt * xm(i, k);
      });

  // Upper (ascending grid) or lower (descending grid) boundary calculation
  Kokkos::parallel_for(
      "windm_edsclrm_rhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        const int i = cell;
        // RHS forcings
        // RHS time tendency
        rhs(i, (nzt - 1)) = real_t(0.5) * (-lhs_diff(2, i, (nzt - 1)) * xm(i, (nzt - 2)) -
                                           lhs_diff(1, i, (nzt - 1)) * xm(i, (nzt - 1))) +
                            xm_tndcy(i, (nzt - 1)) + invrs_dt * xm(i, (nzt - 1));
      });

  if (stats.l_sample) {
    // Native OpenACC host update: the stats API stages its own inputs.
  }

  if (stats.l_sample && var_on_stats_list(stats, name_ta)) {
    Kokkos::parallel_for(
        "windm_edsclrm_rhs::array4", (stats_tmp.extent(0) * stats_tmp.extent(1)),
        KOKKOS_LAMBDA(int cell) {
          const int slice_0           = ((cell % stats_tmp.extent(0)));
          const int slice_1           = ((cell / stats_tmp.extent(0)));
          stats_tmp(slice_0, slice_1) = real_t(0.0);
        });

    Kokkos::parallel_for(
        "windm_edsclrm_rhs::loop5", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = cell;
          int k;
          // Lower boundary
          stats_tmp(i, 0) =
              real_t(0.5) * (lhs_diff(1, i, 0) * xm(i, 0) + lhs_diff(0, i, 0) * xm(i, 1));

          for (k = 1; k < (nzt - 1); k += 1) {
            stats_tmp(i, k) =
                real_t(0.5) *
                (lhs_diff((gr_grid_dir_indx + 1), i, k) * xm(i, (-gr_grid_dir_indx + k)) +
                 lhs_diff(1, i, k) * xm(i, k) +
                 lhs_diff((-gr_grid_dir_indx + 1), i, k) * xm(i, (gr_grid_dir_indx + k)));
          }

          // Upper boundary
          stats_tmp(i, (nzt - 1)) = real_t(0.5) * (lhs_diff(2, i, (nzt - 1)) * xm(i, (nzt - 2)) +
                                                   lhs_diff(1, i, (nzt - 1)) * xm(i, (nzt - 1)));
        });

    // Statistics:  explicit contributions for um or vm.
    // xm term ta has both implicit and explicit components; call
    // stats_begin_budget for the explicit part of turbulent advection.
    stats_begin_budget(name_ta, stats_tmp, // In
                       stats);             // InOut
  }

  if (!l_imp_sfc_momentum_flux) {

    // RHS generalized surface flux.
    Kokkos::parallel_for(
        "windm_edsclrm_rhs::loop6", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i        = cell;
          rhs(i, gr_k_lb_zt) = rhs(i, gr_k_lb_zt) + gr_grid_dir * invrs_rho_ds_zt(i, gr_k_lb_zt) *
                                                        gr_invrs_dzt(i, gr_k_lb_zt) *
                                                        rho_ds_zm(i, gr_k_lb_zm) * xpwp_sfc(i);
        });

    if (stats.l_sample) {
      // Native OpenACC host update: the stats API stages its own inputs.
    }
    if (stats.l_sample && var_on_stats_list(stats, name_ta)) {
      Kokkos::parallel_for(
          "windm_edsclrm_rhs::array7", (stats_tmp.extent(0) * stats_tmp.extent(1)),
          KOKKOS_LAMBDA(int cell) {
            const int slice_0           = ((cell % stats_tmp.extent(0)));
            const int slice_1           = ((cell / stats_tmp.extent(0)));
            stats_tmp(slice_0, slice_1) = real_t(0.0);
          });
      Kokkos::parallel_for(
          "windm_edsclrm_rhs::loop8", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i              = cell;
            stats_tmp(i, gr_k_lb_zt) = gr_grid_dir * invrs_rho_ds_zt(i, gr_k_lb_zt) *
                                       gr_invrs_dzt(i, gr_k_lb_zt) * rho_ds_zm(i, gr_k_lb_zm) *
                                       xpwp_sfc(i);
          });
      stats_update_budget(name_ta, stats_tmp, // In
                          stats);             // InOut
    }

    // l_imp_sfc_momentum_flux
  }

  return;
}

} // namespace
} // namespace clubb_kokkos
