// Case forcing and surface boundary coordination; see prescribe_forcings.hpp.
#include "scratch_workspace.hpp"
#include "prescribe_forcings.hpp"
#include "arm.hpp"
#include "arm_0003.hpp"
#include "arm_3year.hpp"
#include "astex_a209.hpp"
#include "atex.hpp"
#include "atex_long.hpp"
#include "arm_97.hpp"
#include "bomex.hpp"
#include "clex9_nov02.hpp"
#include "clex9_oct14.hpp"
#include "cloud_feedback.hpp"
#include "cobra.hpp"
#include "dycoms2_rf01.hpp"
#include "dycoms2_rf02.hpp"
#include "fire.hpp"
#include "gabls2.hpp"
#include "gabls3.hpp"
#include "gabls3_night.hpp"
#include "jun25.hpp"
#include "lba.hpp"
#include "mpace_a.hpp"
#include "mpace_b.hpp"
#include "nov11.hpp"
#include "rico.hpp"
#include "neutral_case.hpp"
#include "ekman.hpp"
#include "twp_ice.hpp"
#include "wangara.hpp"
#include "sfc_flux.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

static void read_surface_var_for_bc(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                                    Real2D um, Real2D vm, Real2D rtm, Real2D thlm, Real2D rho_zm,
                                    Real2D exner, Real1D p_sfc, bool l_modify_bc_for_cnvg_test,
                                    Real1D z_bot, Real1D um_bot, Real1D vm_bot, Real1D rtm_bot,
                                    Real1D thlm_bot, Real1D rho_bot, Real1D exner_bot);

//----------------------------------------------------------------------
void prescribe_forcings(ScratchWorkspace& workspace,                                 // InOut
                        const grid& gr, int nzm, int nzt, int ngrdcol,               // In
                        int sclr_dim, int edsclr_dim, const sclr_idx_type& sclr_idx, // In
                        bool l_do_expldiff_rtm_thlm,                                 // In
                        const std::string& runtype, int sfctype,                     // In
                        double time_current, double time_initial, real_t dt,         // In
                        bool l_readiopdata,                                          // In
                        Real2D um, Real2D vm, Real2D thlm,                           // In
                        Real2D p_in_Pa, Real2D exner, Real2D rho, Real2D rho_zm, Real2D thvm, // In
                        Real1D veg_T_in_K,                                                    // In
                        bool l_modify_bc_for_cnvg_test, int saturation_formula,               // In
                        stats_type& stats,                                                 // InOut
                        bool l_add_dycore_grid, int grid_remap_method,                     // In
                        int total_idx_rho_lin_spline, Real2D rho_lin_spline_vals,          // In
                        Real2D rho_lin_spline_levels, const grid& gr_dycore,               // In
                        Real2D rtm, Real2D wm_zm, Real2D wm_zt, Real2D ug, Real2D vg,      // InOut
                        Real2D um_ref, Real2D vm_ref,                                      // InOut
                        Real2D thlm_forcing, Real2D rtm_forcing, Real2D um_forcing,        // InOut
                        Real2D vm_forcing, Real2D wprtp_forcing, Real2D wpthlp_forcing,    // InOut
                        Real2D rtp2_forcing, Real2D thlp2_forcing, Real2D rtpthlp_forcing, // InOut
                        Real3D wpsclrp, Real3D sclrm_forcing, Real3D edsclrm_forcing,      // InOut
                        Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D upwp_sfc,
                        Real1D vpwp_sfc,                                                    // InOut
                        Real1D T_sfc, Real1D p_sfc, real_t& sens_ht, real_t& latent_ht,     // InOut
                        Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc, err_info_type& err_info,  // InOut
                        const time_dependent_input_state& input, const mpace_a_state& mpace // In
) {
  // Description:
  // Calculate case tendencies and surface variables in Native stage order.
  // The appended contexts provide case-owned data formerly accessed globally;
  // statistics use the same case-owned C++ context as the driver.
  // Match Native CHARACTER(len=50) and its blank-padded SELECT CASE operands.
  // Only trailing blanks are removed; leading blanks remain significant.
  auto workspace_scope           = workspace.scope();
  std::string case_name = runtype.substr(0, 50);
  const auto last_case_character = case_name.find_last_not_of(' ');
  case_name.erase(last_case_character == std::string::npos ? 0 : (last_case_character + 1));
  Real1D ustar                 = workspace_scope.view<Real1D>("prescribe_forcings::ustar", ngrdcol);
  Real1D ubar                  = workspace_scope.view<Real1D>("prescribe_forcings::ubar", ngrdcol);
  bool l_compute_momentum_flux = false;
  bool l_set_sclr_sfc_rtm_thlm = false;
  bool l_fixed_flux            = false;
  Real1D um_bot    = workspace_scope.view<Real1D>("prescribe_forcings::um_bot", ngrdcol);
  Real1D vm_bot    = workspace_scope.view<Real1D>("prescribe_forcings::vm_bot", ngrdcol);
  Real1D rtm_bot   = workspace_scope.view<Real1D>("prescribe_forcings::rtm_bot", ngrdcol);
  Real1D thlm_bot  = workspace_scope.view<Real1D>("prescribe_forcings::thlm_bot", ngrdcol);
  Real1D rho_bot   = workspace_scope.view<Real1D>("prescribe_forcings::rho_bot", ngrdcol);
  Real1D exner_bot = workspace_scope.view<Real1D>("prescribe_forcings::exner_bot", ngrdcol);
  Real1D z_bot     = workspace_scope.view<Real1D>("prescribe_forcings::z_bot", ngrdcol);
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  //-----------------------------------------------------------------------
  // FIND ALL DIAGNOSTIC VARIABLES
  //-----------------------------------------------------------------------
  // Set vertical velocity and compute large-scale forcing. Reset the forcing
  // arrays each iteration, as moved here from the individual case routines.
  Kokkos::parallel_for(
      "prescribe_forcings::reset_thermo", Policy(0, (ngrdcol * nzt)), KOKKOS_LAMBDA(int index) {
        const int i        = (index % ngrdcol);
        const int k        = (index / ngrdcol);
        rtm_forcing(i, k)  = zero;
        thlm_forcing(i, k) = zero;
      });
  Kokkos::parallel_for(
      "prescribe_forcings::reset_momentum", Policy(0, (ngrdcol * nzm)), KOKKOS_LAMBDA(int index) {
        const int i           = (index % ngrdcol);
        const int k           = (index / ngrdcol);
        wprtp_forcing(i, k)   = zero;
        wpthlp_forcing(i, k)  = zero;
        rtp2_forcing(i, k)    = zero;
        thlp2_forcing(i, k)   = zero;
        rtpthlp_forcing(i, k) = zero;
      });
  if (input.l_t_dependent && !input.l_ignore_forcings) {
    // Time-dependent input includes cloud-feedback, GABLS3/night, ARM variants,
    // TWP-ICE, ASTEX and COBRA. Dycore input is remapped onto each physics grid.
    if (l_add_dycore_grid) {
      apply_time_dependent_forcings_from_dycore(
          workspace, ngrdcol, gr.nzm, gr.nzt, sclr_dim, edsclr_dim, sclr_idx,
          l_do_expldiff_rtm_thlm, gr, gr_dycore, time_current, rtm, rho, exner, grid_remap_method,
          total_idx_rho_lin_spline, rho_lin_spline_vals, rho_lin_spline_levels, p_sfc, thlm_forcing,
          rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing, wm_zt, wm_zm, ug, vg, sclrm_forcing,
          edsclrm_forcing, input);
    } else {
      apply_time_dependent_forcings(workspace, ngrdcol, gr.nzm, gr.nzt, sclr_dim, edsclr_dim,
                                    sclr_idx, l_do_expldiff_rtm_thlm, gr, time_current, rtm, rho,
                                    exner, thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing,
                                    vm_forcing, wm_zt, wm_zm, ug, vg, sclrm_forcing,
                                    edsclrm_forcing, input);
    }
    // V.Larson's top-level zero forcing reduces the sponge damping and its
    // associated cloud-feedback sawtooth noise. Retain this for all time-dependent cases.
    Kokkos::parallel_for(
        "prescribe_forcings::top_forcing", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          rtm_forcing(i, (nzt - 1))  = zero;
          thlm_forcing(i, (nzt - 1)) = zero;
        });
  } else {
    // Legacy method of setting case forcing. The disabled ASTEX tendency
    // branch in Native remains inactive; ASTEX forcing comes from its input.
    if (case_name == "atex") {
      // Native passes p_sfc(1) by Fortran sequence association to an array.
      // C++ passes the complete owning View for that same contiguous sequence.
      atex_tndcy(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time_current, time_initial,
                 rtm, l_add_dycore_grid, grid_remap_method, gr_dycore, total_idx_rho_lin_spline,
                 rho_lin_spline_vals, rho_lin_spline_levels, p_sfc, err_info, wm_zt, wm_zm,
                 thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
      if (clubb_at_least_debug_level_api(0)) {
        const auto codes = to_host_space(workspace, err_info.err_code);
        for (int i = 0; i < ngrdcol; ++i) {
          if (codes(i) == clubb_fatal_error) {
            std::cerr << err_info.err_header_global << '\n';
            std::cerr << "Calling atex_tndcy in prescribe_forcings\n";
            return;
          }
        }
      }
    } else if (case_name == "atex_long") {
      atex_long_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time_current,
                       wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "bomex" || case_name == "bomex_5day_4scam") {
      bomex_tndcy(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, rtm, l_readiopdata,
                  thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "dycoms2_rf01") {
      dycoms2_rf01_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx,
                          thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "dycoms2_rf02") {
      dycoms2_rf02_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, wm_zt, wm_zm,
                          thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "fire" || case_name == "generic" || case_name == "coriolis_test" ||
               case_name == "ekman") {
      // Analytic radiation is computed elsewhere.
      Kokkos::parallel_for(
          "prescribe_forcings::legacy_zero", Policy(0, (gr.nzt * ngrdcol)),
          KOKKOS_LAMBDA(int index) {
            const int i        = (index % ngrdcol);
            const int k        = (index / ngrdcol);
            thlm_forcing(i, k) = real_t(0.0);
            rtm_forcing(i, k)  = real_t(0.0);
          });
    } else if (case_name == "gabls2") {
      gabls2_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time_current, time_initial,
                     wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "lba") {
      lba_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr,
                 thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "mpace_a") {
      mpace_a_tndcy(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, time_current, p_in_Pa,
                    wm_zt, wm_zm, thlm_forcing, rtm_forcing, um_ref, vm_ref, sclrm_forcing,
                    edsclrm_forcing, mpace);
    } else if (case_name == "mpace_b") {
      mpace_b_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, p_in_Pa, thvm,
                      wm_zt, wm_zm, thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "rico") {
      rico_tndcy(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, rtm, exner, thlm_forcing,
                 rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else if (case_name == "neutral") {
      Kokkos::parallel_for(
          "prescribe_forcings::neutral", Policy(0, (gr.nzt * ngrdcol)), KOKKOS_LAMBDA(int index) {
            const int i        = (index % ngrdcol);
            const int k        = (index / ngrdcol);
            thlm_forcing(i, k) = real_t(0.0);
            rtm_forcing(i, k)  = real_t(0.0);
          });
      if (sclr_dim > 0)
        Kokkos::deep_copy(sclrm_forcing, real_t(0.0));
      if (edsclr_dim > 0)
        Kokkos::deep_copy(edsclrm_forcing, real_t(0.0));
    } else if (case_name == "wangara") {
      wangara_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx, wm_zt, wm_zm,
                      thlm_forcing, rtm_forcing, sclrm_forcing, edsclrm_forcing);
    } else {
      std::cerr << "prescribe_forcings: Don't know how to handle LS forcing for runtype: "
                << runtype << '\n';
      throw std::runtime_error("");
    }
  }

  //----------------------------------------------------------------
  // Compute Surface Fluxes
  //----------------------------------------------------------------
  // Derive the bottom-level physical input for surface boundary conditions.
  // The convergence-test option instead samples physical quantities at 25 m.
  read_surface_var_for_bc(workspace, gr, ngrdcol, um, vm, rtm, thlm, rho_zm, exner, p_sfc,
                          l_modify_bc_for_cnvg_test, z_bot, um_bot, vm_bot, rtm_bot, thlm_bot,
                          rho_bot, exner_bot);
  // Boundary conditions for second-order moments.
  compute_ubar(ngrdcol, um_bot, vm_bot, ubar);

  if (case_name == "rico") {
    l_set_sclr_sfc_rtm_thlm = true;
    rico_sfclyr(workspace, ngrdcol, time_current, um_bot, vm_bot, thlm_bot, rtm_bot, z_bot, p_sfc,
                exner_bot, saturation_formula, upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar,
                T_sfc, input);
  } else if (case_name == "gabls3") {
    l_compute_momentum_flux = true;
    gabls3_sfclyr(workspace, ngrdcol, ubar, veg_T_in_K, thlm_bot, rtm_bot, z_bot, exner_bot,
                  wpthlp_sfc, wprtp_sfc, ustar);
  } else if (case_name == "gabls3_night") {
    gabls3_night_sfclyr(workspace, ngrdcol, time_current, um_bot, vm_bot, thlm_bot, rtm_bot, z_bot,
                        upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar, input);
  } else if (case_name == "jun25_altocu") {
    // No surface momentum or heat fluxes for the June 25 altocumulus case.
    Kokkos::deep_copy(ustar, real_t(0.0));
    jun25_altocu_read_t_dependent(workspace, time_current, sens_ht, latent_ht, input);
  } else if (case_name == "cobra") {
    l_compute_momentum_flux = true;
    cobra_sfclyr(workspace, ngrdcol, sclr_dim, edsclr_dim, sclr_idx, l_do_expldiff_rtm_thlm,
                 time_current, z_bot, rho_bot, thlm_bot, ubar, wpthlp_sfc, wprtp_sfc, ustar,
                 wpsclrp_sfc, wpedsclrp_sfc, T_sfc, input);
  } else if (case_name == "clex9_nov02") {
    // No surface momentum or heat fluxes for CLEX9 November 2 altocumulus.
    Kokkos::deep_copy(ustar, real_t(0.0));
    clex9_nov02_read_t_dependent(workspace, time_current, sens_ht, latent_ht, input);
  } else if (case_name == "clex9_oct14") {
    // No surface momentum or heat fluxes for CLEX9 October 14 altocumulus.
    Kokkos::deep_copy(ustar, real_t(0.0));
    clex9_oct14_read_t_dependent(workspace, time_current, sens_ht, latent_ht, input);
  } else if (case_name == "astex_a209") {
    l_compute_momentum_flux = true;
    astex_a209_sfclyr(workspace, ngrdcol, time_current, ubar, rtm_bot, thlm_bot, z_bot, exner_bot,
                      p_sfc, saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, T_sfc, input);
  } else if (case_name == "nov11_altocu") {
    // No surface momentum or heat fluxes. The one-time water adjustment from
    // the obsolete tendency routine must still run each timestep here(E. Hoppe, 2011).
    Kokkos::deep_copy(ustar, real_t(0.0));
    nov11_altocu_rtm_adjust(ngrdcol, gr, time_current, time_initial, dt, rtm);
    nov11_altocu_read_t_dependent(workspace, time_current, sens_ht, latent_ht, input);
  } else if (case_name == "fire" || case_name == "generic") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    l_fixed_flux            = true;
    fire_sfclyr(workspace, ngrdcol, time_current, ubar, p_sfc, thlm_bot, rtm_bot, exner_bot,
                saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, T_sfc, input);
  } else if (case_name == "cloud_feedback_s6" || case_name == "cloud_feedback_s6_p2k" ||
             case_name == "cloud_feedback_s11" || case_name == "cloud_feedback_s11_p2k" ||
             case_name == "cloud_feedback_s12" || case_name == "cloud_feedback_s12_p2k" ||
             case_name == "cgils_s6" || case_name == "cgils_s6_p2k" || case_name == "cgils_s11" ||
             case_name == "cgils_s11_p2k" || case_name == "cgils_s12" || case_name == "cgils_s12_p2k") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    l_fixed_flux            = true;
    cloud_feedback_sfclyr(workspace, ngrdcol, time_current, sfctype, thlm_bot, rtm_bot, z_bot, ubar,
                          p_sfc, T_sfc, saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, input);
  } else if (case_name == "arm") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    arm_sfclyr(workspace, ngrdcol, time_current, z_bot, thlm_bot, ubar, wpthlp_sfc, wprtp_sfc,
               ustar, input);
  } else if (case_name == "arm_0003") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    arm_0003_sfclyr(workspace, ngrdcol, time_current, z_bot, rho_bot, thlm_bot, ubar, wpthlp_sfc,
                    wprtp_sfc, ustar, input);
  } else if (case_name == "arm_3year") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    arm_3year_sfclyr(workspace, ngrdcol, time_current, z_bot, rho_bot, thlm_bot, ubar, wpthlp_sfc,
                     wprtp_sfc, ustar, input);
  } else if (case_name == "arm_97" || case_name == "mc3e") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    arm_97_sfclyr(workspace, ngrdcol, time_current, z_bot, rho_bot, thlm_bot, ubar, wpthlp_sfc,
                  wprtp_sfc, ustar, input);
  } else if (case_name == "atex") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    atex_sfclyr(workspace, ngrdcol, time_current, ubar, thlm_bot, rtm_bot, exner_bot, wpthlp_sfc,
                wprtp_sfc, ustar, T_sfc, input);
  } else if (case_name == "atex_long") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    atex_long_sfclyr(workspace, ngrdcol, time_current, ubar, thlm_bot, rtm_bot, exner_bot, rho_bot,
                     wpthlp_sfc, wprtp_sfc, ustar, T_sfc, input);
  } else if (case_name == "bomex" || case_name == "bomex_5day_4scam") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    bomex_sfclyr(workspace, ngrdcol, time_current, rtm_bot, wpthlp_sfc, wprtp_sfc, ustar, input);
  } else if (case_name == "dycoms2_rf01") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    dycoms2_rf01_sfclyr(workspace, ngrdcol, time_current, sfctype, p_sfc, exner_bot, ubar, thlm_bot,
                        rtm_bot, rho_bot, saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, T_sfc,
                        input);
  } else if (case_name == "dycoms2_rf02") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    dycoms2_rf02_sfclyr(workspace, ngrdcol, time_current, wpthlp_sfc, wprtp_sfc, ustar, input);
  } else if (case_name == "gabls2") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    gabls2_sfclyr(workspace, ngrdcol, time_current, time_initial, z_bot, p_sfc, ubar, thlm_bot,
                  rtm_bot, exner_bot, saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, T_sfc);
  } else if (case_name == "lba") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    lba_sfclyr(ngrdcol, time_current, time_initial, z_bot, rho_bot, thlm_bot, ubar,
                 wpthlp_sfc, wprtp_sfc, ustar);
  } else if (case_name == "mpace_a") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    mpace_a_sfclyr(workspace, ngrdcol, time_current, rho_bot, wpthlp_sfc, wprtp_sfc, ustar, mpace);
  } else if (case_name == "mpace_b") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    mpace_b_sfclyr(workspace, ngrdcol, time_current, rho_bot, wpthlp_sfc, wprtp_sfc, ustar, input);
  } else if (case_name == "neutral") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    neutral_case_sfclyr(ngrdcol, time_current, um_bot, vm_bot, ubar, upwp_sfc, vpwp_sfc,
                          wpthlp_sfc, wprtp_sfc, ustar);
  } else if (case_name == "ekman") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    ekman_sfclyr(ngrdcol, um_bot, vm_bot, ubar, upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar);
  } else if (case_name == "twp_ice") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    twp_ice_sfclyr(workspace, ngrdcol, time_current, z_bot, exner_bot, thlm_bot, ubar, rtm_bot,
                   p_sfc, saturation_formula, wpthlp_sfc, wprtp_sfc, ustar, T_sfc, input);
  } else if (case_name == "wangara") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = true;
    wangara_sfclyr(ngrdcol, time_current, wpthlp_sfc, wprtp_sfc, ustar);
  } else if (case_name == "coriolis_test") {
    l_compute_momentum_flux = true;
    l_set_sclr_sfc_rtm_thlm = false;
    l_fixed_flux            = true;
    Kokkos::deep_copy(ustar, real_t(0.0));
  } else {
    std::cerr << "Invalid value of runtype = " << runtype << '\n';
    throw std::runtime_error("");
  }

  // Shared final operations preserve the case-specific flags above.
  if (l_compute_momentum_flux)
    compute_momentum_flux(ngrdcol, um_bot, vm_bot, ubar, ustar, upwp_sfc, vpwp_sfc);
  if (l_set_sclr_sfc_rtm_thlm)
    set_sclr_sfc_rtm_thlm(ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                           wpthlp_sfc, wprtp_sfc, wpsclrp_sfc, wpedsclrp_sfc);

  // Fixed-flux surface type overrides the physical fluxes only for flagged cases.
  if (sfctype == 0 && l_fixed_flux) {
    Kokkos::parallel_for(
        "prescribe_forcings::fixed_flux", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          wpthlp_sfc(i) = sens_ht;
          wprtp_sfc(i)  = latent_ht;
        });
    if (sclr_idx.iisclr_thl > -1) {
      Kokkos::parallel_for(
          "prescribe_forcings::fixed_theta_scalar", Policy(0, (ngrdcol * nzm)),
          KOKKOS_LAMBDA(int index) {
            wpsclrp((index % ngrdcol), (index / ngrdcol), sclr_idx.iisclr_thl) = sens_ht;
          });
    }
    if (sclr_idx.iisclr_rt > -1) {
      Kokkos::parallel_for(
          "prescribe_forcings::fixed_water_scalar", Policy(0, (ngrdcol * nzm)),
          KOKKOS_LAMBDA(int index) {
            wpsclrp((index % ngrdcol), (index / ngrdcol), sclr_idx.iisclr_rt) = latent_ht;
          });
    }
  }

  // Save the surface fluxes through the retained statistics backend. Its public
  // array boundary owns host staging; no statistics state is copied.
  if (stats.l_sample) {
    Real1D sh = workspace_scope.view<Real1D>("prescribe_forcings::sh", ngrdcol);
    Real1D lh = workspace_scope.view<Real1D>("prescribe_forcings::lh", ngrdcol);
    Kokkos::parallel_for(
        "prescribe_forcings::sensible_stats", Policy(0, ngrdcol),
        KOKKOS_LAMBDA(int i) { sh(i) = wpthlp_sfc(i) * rho_zm(i, 0) * Cp; });
    stats_update("sh", sh, stats);
    Kokkos::parallel_for(
        "prescribe_forcings::latent_stats", Policy(0, ngrdcol),
        KOKKOS_LAMBDA(int i) { lh(i) = wprtp_sfc(i) * rho_zm(i, 0) * Lv; });
    stats_update("lh", lh, stats);
    stats_update("wpthlp_sfc", wpthlp_sfc, stats);
    stats_update("wprtp_sfc", wprtp_sfc, stats);
    stats_update("upwp_sfc", upwp_sfc, stats);
    stats_update("vpwp_sfc", vpwp_sfc, stats);
    stats_update("ustar", ustar, stats);
    stats_update("T_sfc", T_sfc, stats);
  }
}

//-------------------------------------------------------------------------------
static void read_surface_var_for_bc(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                                    Real2D um, Real2D vm, Real2D rtm, Real2D thlm, Real2D rho_zm,
                                    Real2D exner, Real1D p_sfc, bool l_modify_bc_for_cnvg_test,
                                    Real1D z_bot, Real1D um_bot, Real1D vm_bot, Real1D rtm_bot,
                                    Real1D thlm_bot, Real1D rho_bot, Real1D exner_bot) {
  // Description:
  // Derive bottom-level physical quantities for surface boundary conditions.
  // The default uses the lowest model values; the convergence-test option
  // interpolates to a fixed 25 m height to remove vertical-grid dependence from
  // the surface-flux calculation. Author: Shixuan Zhang (PNNL).
  auto workspace_scope                 = workspace.scope();
  constexpr int constant_height_option = 2; //1: nearest level; 2: fixed-height interpolation.
  const int nzt                        = gr.nzt;
  const int nzm                        = gr.nzm;
  const auto zt                        = gr.zt;
  const auto zm                        = gr.zm;
  Real2D um_zm = (l_modify_bc_for_cnvg_test && constant_height_option != 1)
                     ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                    "read_surface_var_for_bc::um_zm", ngrdcol, nzm)
                     : Real2D{};
  Real2D vm_zm = (l_modify_bc_for_cnvg_test && constant_height_option != 1)
                     ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                    "read_surface_var_for_bc::vm_zm", ngrdcol, nzm)
                     : Real2D{};
  Real2D exner_zm =
      (l_modify_bc_for_cnvg_test && constant_height_option != 1)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "read_surface_var_for_bc::exner_zm", ngrdcol, nzm)
          : Real2D{};
  Real2D rtm_zm =
      (l_modify_bc_for_cnvg_test && constant_height_option != 1)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "read_surface_var_for_bc::rtm_zm", ngrdcol, nzm)
          : Real2D{};
  Real2D thlm_zm =
      (l_modify_bc_for_cnvg_test && constant_height_option != 1)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "read_surface_var_for_bc::thlm_zm", ngrdcol, nzm)
          : Real2D{};
  using Levels = Kokkos::View<int*, Kokkos::LayoutLeft>;
  Levels k_min = l_modify_bc_for_cnvg_test
                     ? workspace_scope.view<Levels>(Kokkos::WithoutInitializing,
                                                    "read_surface_var_for_bc::k_min", ngrdcol)
                     : Levels{};
  using Policy = Kokkos::RangePolicy<typename Real2D::execution_space>;

  if (!l_modify_bc_for_cnvg_test) {
    // Default CLUBB-SCM setup. Surface Exner comes from surface pressure.
    Kokkos::parallel_for(
        "read_surface_var_for_bc::default", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          z_bot(i)     = zt(i, 0);
          um_bot(i)    = um(i, 0);
          vm_bot(i)    = vm(i, 0);
          rtm_bot(i)   = rtm(i, 0);
          thlm_bot(i)  = thlm(i, 0);
          rho_bot(i)   = rho_zm(i, 0);
          exner_bot(i) = Kokkos::pow(p_sfc(i) / p0, kappa);
        });
  } else {
    // Modified option: find the closest model level to the fixed 25 m height.
    Kokkos::parallel_for(
        "read_surface_var_for_bc::nearest", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
          z_bot(i)       = real_t(25.0);
          real_t min_val = Kokkos::abs(zt(i, 0) - z_bot(i));
          k_min(i)       = 0;
          for (int k = 1; k < nzt; ++k) {
            if (Kokkos::abs(zt(i, k) - z_bot(i)) < min_val) {
              min_val  = Kokkos::abs(zt(i, k) - z_bot(i));
              k_min(i) = k;
            }
          }
        });
    if (clubb_at_least_debug_level_api(1)) {
      const auto nearest = to_host_space(workspace, k_min);
      for (int i = 0; i < ngrdcol; ++i) {
        if (nearest(i) < 0 || nearest(i) >= nzt) {
          std::cerr << "Sanity check failed! constant model height is not properly set\n";
          std::cerr << "in get_fixed_height_values at i = " << i + 1 << '\n';
        }
      }
    }
    if (constant_height_option == 1) {
      // Source option 1: select the nearest level without interpolation.
      Kokkos::parallel_for(
          "read_surface_var_for_bc::nearest_values", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
            um_bot(i)    = um(i, k_min(i));
            vm_bot(i)    = vm(i, k_min(i));
            rtm_bot(i)   = rtm(i, k_min(i));
            thlm_bot(i)  = thlm(i, k_min(i));
            rho_bot(i)   = rho_zm(i, k_min(i));
            exner_bot(i) = exner(i, k_min(i));
          });
    } else {
      // Source option 2: interpolate all five fields to momentum levels first,
      // then apply the monotone cubic at 25 m. Keep the original call ordering.
      zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, um, um_zm);
      zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, vm, vm_zm);
      zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, thlm, thlm_zm);
      zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, rtm, rtm_zm);
      zt2zm_api(workspace, nzm, nzt, ngrdcol, gr, exner, exner_zm);
      Kokkos::parallel_for(
          "read_surface_var_for_bc::interpolate", Policy(0, ngrdcol), KOKKOS_LAMBDA(int i) {
            exner_zm(i, 0) = Kokkos::pow(p_sfc(i) / p0, kappa);
            int km1, k00, kp1, kp2; // Zero-based interpolation stencil indices.
            if (k_min(i) == 0) {
              km1 = 0;
              k00 = 0;
              kp1 = 1;
              kp2 = 2;
            } else {
              km1 = k_min(i) - 1;
              k00 = k_min(i);
              kp1 = k_min(i) + 1;
              kp2 = k_min(i) + 2;
            }
            rho_bot(i)   = rho_zm(i, k_min(i));
            um_bot(i)    = mono_cubic_interp(z_bot(i), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00),
                                             zm(i, kp1), zm(i, kp2), um_zm(i, km1), um_zm(i, k00),
                                             um_zm(i, kp1), um_zm(i, kp2));
            vm_bot(i)    = mono_cubic_interp(z_bot(i), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00),
                                             zm(i, kp1), zm(i, kp2), vm_zm(i, km1), vm_zm(i, k00),
                                             vm_zm(i, kp1), vm_zm(i, kp2));
            exner_bot(i) = mono_cubic_interp(z_bot(i), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00),
                                             zm(i, kp1), zm(i, kp2), exner_zm(i, km1),
                                             exner_zm(i, k00), exner_zm(i, kp1), exner_zm(i, kp2));
            thlm_bot(i)  = mono_cubic_interp(z_bot(i), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00),
                                             zm(i, kp1), zm(i, kp2), thlm_zm(i, km1),
                                             thlm_zm(i, k00), thlm_zm(i, kp1), thlm_zm(i, kp2));
            rtm_bot(i)   = mono_cubic_interp(z_bot(i), km1, k00, kp1, kp2, zm(i, km1), zm(i, k00),
                                             zm(i, kp1), zm(i, kp2), rtm_zm(i, km1), rtm_zm(i, k00),
                                             rtm_zm(i, kp1), rtm_zm(i, kp2));
          });
    }
  }
}

} // namespace clubb_kokkos
