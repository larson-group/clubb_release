#include "scratch_workspace.hpp"
#include "CLUBB_core/stats_clubb_utilities.hpp"
// Microphysics calculations translated from the same-named Native source.
//
// Public coordinators take caller-owned column-first LayoutLeft Views and retain
// Native call, statistics and recurrence order. Scalar helpers execute on the
// Kokkos backend; scratch allocations live only for the current call.
// Local storage coordinates, sedimentation level arguments and cloud-top results
// are zero-based. Shared hydrometeor metadata IDs retain their Native contract.
#pragma once
#include "Microphys/microphys_diagnostics.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/index_mapping.hpp"
#include "CLUBB_core/T_in_K_module.hpp"
#include "CLUBB_core/saturation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/fill_holes.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/matrix_solver_wrapper.hpp"
#include "CLUBB_core/diffusion.hpp"
#include "CLUBB_core/mean_adv.hpp"
#include "Microphys/parameters_microphys.hpp"
namespace clubb_kokkos {
using Int1D = Kokkos::View<int*, Kokkos::LayoutLeft>;
using Int3D = Kokkos::View<int***, Kokkos::LayoutLeft>;
inline void advance_microphys(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, double time_current,
    int hydromet_dim, const hm_metadata_type& hm_metadata, Real2D wm_zt, Real2D wp2, Real2D exner,
    Real2D rho, Real2D rho_zm, Real2D rcm, Real2D cloud_frac, Real2D Kh_zm, Real2D Skw_zm,
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real3D hydromet_mc, Real2D Ncm_mc,
    Real2D Lscale, Real3D hydromet_vel_covar_zt_impc, Real3D hydromet_vel_covar_zt_expc,
    Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep, int tridiag_solve_method,
    int fill_holes_type, bool l_upwind_xm_ma, stats_type& stats, Real3D hydromet,
    Real3D hydromet_vel_zt, Real3D hydrometp2, Real3D K_hm, Real2D Ncm, Real2D Nc_in_cloud,
    Real2D rvm_mc, Real2D thlm_mc, err_info_type& err_info, Real3D wphydrometp, Real2D wpNcp,
    const parameters_microphys& microphys);
inline void advance_hydrometeor(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, Real2D wm_zt, Real2D exner, Real2D cloud_frac, Real3D K_hm,
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real3D hydromet_mc,
    Real3D hydromet_vel_covar_zt_impc, Real3D hydromet_vel_covar_zt_expc,
    const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma, int tridiag_solve_method,
    int fill_holes_type, stats_type& stats, Real3D hydromet, Real3D hydromet_vel_zt,
    Real3D hydrometp2, Real2D rvm_mc, Real2D thlm_mc, err_info_type& err_info, Real3D wphydrometp,
    Real3D hydromet_vel, Real3D hydromet_vel_covar, Real3D hydromet_vel_covar_zt,
    const parameters_microphys& microphys);
inline void advance_Ncm(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                        Real2D wm_zt, Real2D cloud_frac, Real2D K_Nc, Real2D rcm, Real2D rho_ds_zm,
                        Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real2D Ncm_mc,
                        const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma,
                        int tridiag_solve_method, stats_type& stats, Real2D Ncm, Real2D Nc_in_cloud,
                        err_info_type& err_info, Real2D wpNcp,
                        const parameters_microphys& microphys);
inline void microphys_solve(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                            const std::string& solve_type, bool l_sed, Real3D lhs_ta, Real3D lhs_ma,
                            Real3D sed_turb_lhs, Real3D sed_diff_lhs, Real2D cloud_frac,
                            int tridiag_solve_method, stats_type& stats, Real3D lhs, Real2D rhs,
                            Real2D hmm, err_info_type& err_info,
                            const parameters_microphys& microphys);
inline void microphys_lhs(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                          const std::string& solve_type, bool l_sed, real_t dt, Real2D K_hm,
                          Real1D nu, Real2D wm_zt, Real2D V_hm, Real2D V_hmt,
                          Real2D Vhmphmp_zt_impc, Real2D rho_ds_zm, Real2D rho_ds_zt,
                          Real2D invrs_rho_ds_zt, bool l_upwind_xm_ma, stats_type& stats,
                          Real3D lhs_ta, Real3D lhs_ma, Real3D sed_turb_lhs, Real3D sed_diff_lhs,
                          Real3D lhs, const parameters_microphys& microphys);
inline void microphys_rhs(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                          const std::string& solve_type, real_t dt, bool l_sed, Real2D hmm,
                          Real2D hmm_tndcy, Real2D K_hm, Real1D nu, Real2D cloud_frac,
                          Real2D Vhmphmp_zt_expc, Real2D rho_ds_zm, Real2D rho_ds_zt,
                          Real2D invrs_rho_ds_zt, stats_type& stats, Real2D rhs,
                          const parameters_microphys& microphys);
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
sed_centered_diff_lhs(const grid& gr, real_t V_hmp1, real_t V_hm, real_t rho_ds_zmp1,
                      real_t rho_ds_zm, real_t invrs_rho_ds_zt, real_t invrs_dzt, int level);
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
sed_upwind_diff_lhs(const grid& gr, real_t V_hmt, real_t V_hmtp1, real_t rho_ds_zt,
                    real_t rho_ds_ztp1, real_t invrs_rho_ds_zt, real_t invrs_dzmp1, int level);
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
term_turb_sed_lhs(const grid& gr, real_t Vhmphmp_impcp1, real_t Vhmphmp_impc,
                  real_t Vhmphmp_zt_impcp1, real_t Vhmphmp_zt_impc, real_t rho_ds_zmp1,
                  real_t rho_ds_zm, real_t rho_ds_ztp1, real_t rho_ds_zt, real_t invrs_dzt,
                  real_t invrs_dzmp1, real_t invrs_rho_ds_zt, int level, bool l_upwind_diff_sed);
KOKKOS_INLINE_FUNCTION real_t term_turb_sed_rhs(const grid& gr, real_t Vhmphmp_expcp1,
                                                real_t Vhmphmp_expc, real_t Vhmphmp_zt_expcp1,
                                                real_t Vhmphmp_zt_expc, real_t rho_ds_zmp1,
                                                real_t rho_ds_zm, real_t rho_ds_ztp1,
                                                real_t rho_ds_zt, real_t invrs_dzt,
                                                real_t invrs_dzmp1, real_t invrs_rho_ds_zt,
                                                int level, bool l_upwind_diff_sed);
inline Real3D calculate_K_hm(const grid& gr, int ngrdcol, Real2D wp2, Real2D Kh_zm, Real2D Skw_zm,
                             Real2D Lscale, int hydromet_dim, Real1D hydromet_tol, Real3D hydromet,
                             Real3D hydrometp2, Real2D clubb_params, bool l_use_non_local_diff_fac);
inline Int1D get_cloud_top_level(int nzt, int ngrdcol, Real2D rcm, Real3D hydromet,
                                 int hydromet_dim, int iiri);
inline void get_cloud_top_level(int nzt, int ngrdcol, Real2D rcm, Real3D hydromet, int hydromet_dim,
                                int iiri, Int1D cloud_top_level);
inline void write_adv_micro_errors(
    const grid& gr, int ngrdcol, real_t dt, double time_current, int hydromet_dim, Real2D wm_zt,
    Real2D wp2, Real2D exner, Real2D rho, Real2D rho_zm, Real2D rcm, Real2D cloud_frac,
    Real2D Kh_zm, Real2D Skw_zm, Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt,
    Real3D hydromet_mc, Real2D Ncm_mc, Real2D Lscale, Real3D hydromet_vel_covar_zt_impc,
    Real3D hydromet_vel_covar_zt_expc, Real2D clubb_params,
    const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma, Real3D hydromet,
    Real3D hydromet_vel_zt, Real3D hydrometp2, Real3D K_hm, Real2D Ncm, Real2D Nc_in_cloud,
    Real2D rvm_mc, Real2D thlm_mc, Real3D wphydrometp, Real2D wpNcp, const err_info_type& err_info);
inline void advance_microphys(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, double time_current,
    int hydromet_dim, const hm_metadata_type& hm_metadata, Real2D wm_zt, Real2D wp2, Real2D exner,
    Real2D rho, Real2D rho_zm, Real2D rcm, Real2D cloud_frac, Real2D Kh_zm, Real2D Skw_zm,
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real3D hydromet_mc, Real2D Ncm_mc,
    Real2D Lscale, Real3D hydromet_vel_covar_zt_impc, Real3D hydromet_vel_covar_zt_expc,
    Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep, int tridiag_solve_method,
    int fill_holes_type, bool l_upwind_xm_ma, stats_type& stats, Real3D hydromet,
    Real3D hydromet_vel_zt, Real3D hydrometp2, Real3D K_hm, Real2D Ncm, Real2D Nc_in_cloud,
    Real2D rvm_mc, Real2D thlm_mc, err_info_type& err_info, Real3D wphydrometp, Real2D wpNcp,
    const parameters_microphys& microphys) {
  // Description:
  // Advance mean precipitating hydrometeors and mean cloud droplet
  // concentration one model time step, and calculate some statistics.
  // References:
  //---------------------------------------------------------------------------
  //---------------------- Input Variables ----------------------
  // approximation rather than a centered differencing for turbulent or
  // mean advection terms. It affects rtm, thlm, sclrm, um and vm.
  //---------------------- Output Variables ----------------------
  //---------------------- Local Variables ----------------------
  // Turbulent advection for hydrometeors -- down-gradient approximation for
  // covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
  // <w'hm'> = - K_hm * d<hm>/dz; and
  // <w'Nc'> = - K_Nc & d<Nc>/dz;
  // where the coefficients of diffusion, K_hm and K_Nc are variable and depend
  // on multiple factors.
  // eddy-diffusivity applied to
  // hydrometeors
  // prevent turbulent advection of
  // hydrometeors to those levels.
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope                      = workspace.scope();
  const auto l_predict_Nc                   = microphys.l_predict_Nc;
  const auto microphys_scheme               = microphys.microphys_scheme;
  const auto microphys_start_time           = microphys.microphys_start_time;
  const auto hm_metadata_hydromet_list_host = to_host_space(workspace, hm_metadata.hydromet_list);
  const bool l_sample                       = stats.l_sample;

  Real3D hydromet_vel = workspace_scope.view<Real3D>(
      "advance_microphys::hydromet_vel", ngrdcol, gr.nzm,
      hydromet_dim); // Mean hydrometeor sedimentation velocity, <V_xx> [m/s]
  Real3D hydromet_vel_covar =
      workspace_scope.view<Real3D>("advance_microphys::hydromet_vel_covar", ngrdcol, gr.nzm,
                                   hydromet_dim); // Covariance of V_xx & x_x (m-levs)  [units(m/s)]
  Real3D hydromet_vel_covar_zt =
      workspace_scope.view<Real3D>("advance_microphys::hydromet_vel_covar_zt", ngrdcol, gr.nzt,
                                   hydromet_dim); // Covariance of V_xx & x_x (t-levs)  [units(m/s)]
  Real2D K_Nc = workspace_scope.view<Real2D>(
      "advance_microphys::K_Nc", ngrdcol,
      gr.nzm); // Coefficient of diffusion (turb. adv.) for Nc             [m^2/s]
  Real2D rrm_zm = workspace_scope.view<Real2D>(
      "advance_microphys::rrm_zm", ngrdcol,
      gr.nzm); // Rain water mixing ratio on momentum levels               [kg/kg]
  int k;       // Loop indices
  Int1D cloud_top_level = workspace_scope.view<Int1D>(
      "advance_microphys::cloud_top_level", ngrdcol); // Vertical level index of cloud top    [-]
  Real1D precip_rate_sfc =
      workspace_scope.view<Real1D>("advance_microphys::precip_rate_sfc", ngrdcol);
  Real1D rain_flux_sfc = workspace_scope.view<Real1D>("advance_microphys::rain_flux_sfc", ngrdcol);
  Real1D rrm_sfc       = workspace_scope.view<Real1D>("advance_microphys::rrm_sfc", ngrdcol);
  constexpr bool l_use_non_local_diff_fac    = false;
  constexpr bool l_prevent_hm_ta_above_cloud = false;
  int iirr; // Variable(s)
  int iiri;
  int iiNi;
  std::string khm_name;
  const auto gr_nzm = gr.nzm;

  //---------------------- Begin Code ----------------------

  iirr = hm_metadata.iirr;
  iiri = hm_metadata.iiri;
  iiNi = hm_metadata.iiNi;

  // Initialize intent(out) variables -- covariances <w'hm'> (for any
  // hydrometeor, hm) and <w'Nc'>.
  if ((hydromet_dim > 0)) {
    Kokkos::parallel_for(
        "advance_microphys::array1",
        ((wphydrometp.extent_int(0) * wphydrometp.extent_int(1) * wphydrometp.extent_int(2))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % wphydrometp.extent_int(0)));
          const int index_1 = (((cell / wphydrometp.extent_int(0)) % wphydrometp.extent_int(1)));
          const int index_2 = ((cell / (wphydrometp.extent_int(0) * wphydrometp.extent_int(1))));
          wphydrometp(index_0, index_1, index_2) = zero;
        });
  }
  Kokkos::parallel_for(
      "advance_microphys::array2", ((wpNcp.extent_int(0) * wpNcp.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0       = ((cell % wpNcp.extent_int(0)));
        const int index_1       = ((cell / wpNcp.extent_int(0)));
        wpNcp(index_0, index_1) = zero;
      });

  // Return if there is delay between the model start time and start of the
  // microphysics
  if ((time_current < microphys_start_time)) {
    return;
  }

  // Turbulent advection for hydrometeors -- down-gradient approximation for
  // covariances <w'hm'> (for any hydrometeor, hm) and <w'Nc'>:
  // <w'hm'> = - K_hm * d<hm>/dz; and
  // <w'Nc'> = - K_Nc & d<Nc>/dz;
  // where the coefficients of diffusion, K_hm and K_Nc are variable and depend
  // on multiple factors.
  if ((hydromet_dim > 0)) {

    // Prevent the turbulent advection of hydrometeors to altitudes above
    // cloud top.
    if (l_prevent_hm_ta_above_cloud) {

      // Find the vertical level index of cloud top.
      get_cloud_top_level(gr.nzt, ngrdcol, rcm, hydromet, hydromet_dim, hm_metadata.iiri,
                          cloud_top_level);

      // l_prevent_hm_ta_above_cloud
    }

    // Solve for the value of K_hm, the coefficient of diffusion for
    // hydrometeors.
    Kokkos::deep_copy(K_hm, calculate_K_hm(gr, ngrdcol, wp2, Kh_zm, Skw_zm, Lscale, hydromet_dim,
                                           hm_metadata.hydromet_tol, hydromet, hydrometp2,
                                           clubb_params, l_use_non_local_diff_fac));

    Kokkos::parallel_for(
        "advance_microphys::loop4", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < hydromet_dim; ++k) {

            // Prevent the turbulent advection of hydrometeors to altitudes above
            // cloud top.
            if (l_prevent_hm_ta_above_cloud) {

              // Set K_hm to 0 above cloud top.  Since K_hm is a momentum-level
              // variable, and since momentum grid levels are located above their
              // corresponding thermodynamic grid levels, set K_hm to 0 starting
              // at the vertical-level index of cloud_top_level.
              if ((((cloud_top_level(i) > 0) && (k != iiri)) && (k != iiNi))) {
                for (int index_0 = 0; index_0 < (gr_nzm - cloud_top_level(i)); ++index_0) {
                  K_hm(i, (cloud_top_level(i) + index_0), k) = zero;
                }
                // cloud_top_level > 0
              }

              // l_prevent_hm_ta_above_cloud
            }

            // Native loop: k = 1, hydromet_dim
          }
        });

    // hydromet_dim > 0
  }

  if (stats.l_sample) {
    for (k = 0; k < hydromet_dim; ++k) {
      if ((static_cast<int>(microphys_detail::trim(
                                std::string(hm_metadata_hydromet_list_host.data() +
                                                ((hm_metadata_hydromet_list_host.extent(0) * k)),
                                            hm_metadata_hydromet_list_host.extent(0)))
                                .size()) > 0)) {
        khm_name = (std::string("K_hm_") +
                    microphys_detail::trim(
                        std::string(hm_metadata_hydromet_list_host.data() +
                                        ((hm_metadata_hydromet_list_host.extent(0) * k)),
                                    hm_metadata_hydromet_list_host.extent(0))
                            .substr(0, 2)));
        stats_update(khm_name, Kokkos::subview(K_hm, Kokkos::ALL, Kokkos::ALL, k), stats);
      }
    }
  }

  if (l_predict_Nc) {

    // Solve for the value of K_Nc, the coefficient of diffusion for cloud
    // droplet concentration.
    Kokkos::parallel_for(
        "advance_microphys::loop5", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzm; ++k) {
            K_Nc(i, k) = (clubb_params(i, ic_K_hm) * Kh_zm(i, k));
            // Native loop: k = 1, gr%nzm, 1
          }
        });

    // l_predict_Nc
  }

  //------------------------------------------------------------------------
  // Advance predictive mean precipitating hydrometeors one model timestep.
  // Loop over all hydrometeor species, apply mean advection, turbulent
  // advection (diffusion), mean sedimentation, turbulent sedimentation, and
  // microphysics tendency terms.
  //------------------------------------------------------------------------

  if ((hydromet_dim > 0)) {

    advance_hydrometeor(workspace, gr, ngrdcol, dt, hydromet_dim, hm_metadata, wm_zt, exner,
                        cloud_frac, K_hm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, hydromet_mc,
                        hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, nu_vert_res_dep,
                        l_upwind_xm_ma, tridiag_solve_method, fill_holes_type, stats, hydromet,
                        hydromet_vel_zt, hydrometp2, rvm_mc, thlm_mc, err_info, wphydrometp,
                        hydromet_vel, hydromet_vel_covar, hydromet_vel_covar_zt, microphys);

    // hydromet_dim > 0
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << std::string("calling advance_hydrometeor") << "\n";
      write_adv_micro_errors(gr, ngrdcol, dt, time_current, hydromet_dim, wm_zt, wp2, exner, rho,
                             rho_zm, rcm, cloud_frac, Kh_zm, Skw_zm, rho_ds_zm, rho_ds_zt,
                             invrs_rho_ds_zt, hydromet_mc, Ncm_mc, Lscale,
                             hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, clubb_params,
                             nu_vert_res_dep, l_upwind_xm_ma, hydromet, hydromet_vel_zt, hydrometp2,
                             K_hm, Ncm, Nc_in_cloud, rvm_mc, thlm_mc, wphydrometp, wpNcp, err_info);
      return;
      // err_info%err_code == clubb_fatal_error
    }
  }

  //-----------------------------------------------------------------------
  // When mean cloud droplet concentration, Ncm, is predicted, apply
  // sedimentation, advection, and diffusion, and advance Ncm one model
  // timestep.
  //-----------------------------------------------------------------------

  if (l_predict_Nc) {

    // Nc is predicted.
    advance_Ncm(workspace, gr, ngrdcol, dt, wm_zt, cloud_frac, K_Nc, rcm, rho_ds_zm, rho_ds_zt,
                invrs_rho_ds_zt, Ncm_mc, nu_vert_res_dep, l_upwind_xm_ma, tridiag_solve_method,
                stats, Ncm, Nc_in_cloud, err_info, wpNcp, microphys);

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << std::string("in advance_Ncm") << "\n";
        write_adv_micro_errors(
            gr, ngrdcol, dt, time_current, hydromet_dim, wm_zt, wp2, exner, rho, rho_zm, rcm,
            cloud_frac, Kh_zm, Skw_zm, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, hydromet_mc, Ncm_mc,
            Lscale, hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc, clubb_params,
            nu_vert_res_dep, l_upwind_xm_ma, hydromet, hydromet_vel_zt, hydrometp2, K_hm, Ncm,
            Nc_in_cloud, rvm_mc, thlm_mc, wphydrometp, wpNcp, err_info);
        return;
      }
    }

  } else {

    // Nc is prescribed.
    // The in-cloud mean of cloud droplet concentration, Nc_in_cloud, is
    // constant in this scenario.  The overall mean of cloud droplet
    // concentration, Ncm, depends of the in-cloud mean and cloud fraction.

    Kokkos::parallel_for(
        "advance_microphys::array6", ((Ncm.extent_int(0) * Ncm.extent_int(1))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0     = ((cell % Ncm.extent_int(0)));
          const int index_1     = ((cell / Ncm.extent_int(0)));
          Ncm(index_0, index_1) = (Nc_in_cloud(index_0, index_1) * cloud_frac(index_0, index_1));
        });

    // l_predict_Nc
  }

  if (stats.l_sample) {
    stats_update(std::string("Ncm"), Ncm, stats);
    stats_update(std::string("Nc_in_cloud"), Nc_in_cloud, stats);
  }

  if ((iirr >= 0)) {

    // Rainfall rate (mm/day) is defined on thermodynamic levels.  The rainfall
    // rate is given by the equation:
    // Rainfall rate (mm/day) = - V_rr * r_r * ( rho_a / rho_lw )
    //                            * ( 86400 s/day ) * ( 1000 mm/m ).
    // The level average rainfall rate is given by:
    // < Rainfall rate (mm/day) > = - < V_rr * r_r > * ( rho_a / rho_lw )
    //                                * ( 86400 s/day ) * ( 1000 mm/m );
    // which can also be written as:
    // < Rainfall rate (mm/day) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
    //                                * ( rho_a / rho_lw )
    //                                * ( 86400 s/day ) * ( 1000 mm/m ).
    // Rainfall rate is defined as positive.  Since V_rr is always negative,
    // the minus (-) sign is necessary.

    if (stats.l_sample) {
      Real2D array_actual_7 = workspace_scope.view<Real2D>(
          "advance_microphys::array_actual_7", hydromet.extent_int(0), hydromet.extent_int(1));
      Kokkos::parallel_for(
          "advance_microphys::array8", ((hydromet.extent_int(0) * hydromet.extent_int(1))),
          KOKKOS_LAMBDA(int cell) {
            const int index_0 = ((cell % hydromet.extent_int(0)));
            const int index_1 = ((cell / hydromet.extent_int(0)));
            array_actual_7(index_0, index_1) =
                ((((((-((hydromet(index_0, index_1, iirr) *
                         hydromet_vel_zt(index_0, index_1, iirr)) +
                        hydromet_vel_covar_zt(index_0, index_1, iirr)))) > (zero)
                        ? ((-((hydromet(index_0, index_1, iirr) *
                               hydromet_vel_zt(index_0, index_1, iirr)) +
                              hydromet_vel_covar_zt(index_0, index_1, iirr))))
                        : (zero)) *
                   (rho(index_0, index_1) / rho_lw)) *
                  sec_per_day) *
                 mm_per_m);
          });
      stats_update(std::string("precip_rate_zt"), array_actual_7, stats);
    }

    // Precipitation Flux (W/m^2) is defined on momentum levels.  The
    // precipitation flux is given by the equation:
    // Precip. flux (W/m^2) = - V_rr * r_r * rho_a * L_v.
    // The level average precipitation flux is given by:
    // < Precip. flux (W/m^2) > = - < V_rr * r_r > * rho_a * L_v;
    // which can also be written as:
    // < Precip. flux (W/m^2) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
    //                              * rho_a * L_v.
    // It is generally a convention in meteorology to show Precipitation Flux
    // as a positive downward quantity, so the minus (-) sign is necessary.

    if (stats.l_sample) {
      zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr,
                Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, iirr), rrm_zm);
      Real2D array_actual_9 = workspace_scope.view<Real2D>(
          "advance_microphys::array_actual_9", rrm_zm.extent_int(0), rrm_zm.extent_int(1));
      Kokkos::parallel_for(
          "advance_microphys::array10", ((rrm_zm.extent_int(0) * rrm_zm.extent_int(1))),
          KOKKOS_LAMBDA(int cell) {
            const int index_0 = ((cell % rrm_zm.extent_int(0)));
            const int index_1 = ((cell / rrm_zm.extent_int(0)));
            array_actual_9(index_0, index_1) =
                (((((-((rrm_zm(index_0, index_1) * hydromet_vel(index_0, index_1, iirr)) +
                       hydromet_vel_covar(index_0, index_1, iirr)))) > (zero)
                       ? ((-((rrm_zm(index_0, index_1) * hydromet_vel(index_0, index_1, iirr)) +
                             hydromet_vel_covar(index_0, index_1, iirr))))
                       : (zero)) *
                  rho_zm(index_0, index_1)) *
                 Lv);
          });
      stats_update(std::string("Fprec"), array_actual_9, stats);
    }

    // Store values of surface fluxes for statistics
    // See notes above.

    if ((microphys_detail::trim(microphys_scheme) != std::string("morrison"))) {
      if (stats.l_sample) {
        Kokkos::parallel_for(
            "advance_microphys::loop11", ngrdcol, KOKKOS_LAMBDA(int column) {
              const int i = column;
              precip_rate_sfc(i) =
                  ((((((-((hydromet(i, 0, iirr) * hydromet_vel_zt(i, 0, iirr)) +
                          hydromet_vel_covar_zt(i, 0, iirr)))) > (zero)
                          ? ((-((hydromet(i, 0, iirr) * hydromet_vel_zt(i, 0, iirr)) +
                                hydromet_vel_covar_zt(i, 0, iirr))))
                          : (zero)) *
                     (rho(i, 0) / rho_lw)) *
                    sec_per_day) *
                   mm_per_m);
            });
        stats_update(std::string("precip_rate_sfc"), precip_rate_sfc, stats);
      }
      // microphys_scheme /= "morrison"
    }

    if (stats.l_sample) {
      Kokkos::parallel_for(
          "advance_microphys::loop12", ngrdcol, KOKKOS_LAMBDA(int column) {
            const int i      = column;
            rain_flux_sfc(i) = (((((-((rrm_zm(i, 0) * hydromet_vel(i, 0, iirr)) +
                                      hydromet_vel_covar(i, 0, iirr)))) > (zero)
                                      ? ((-((rrm_zm(i, 0) * hydromet_vel(i, 0, iirr)) +
                                            hydromet_vel_covar(i, 0, iirr))))
                                      : (zero)) *
                                 rho_zm(i, 0)) *
                                Lv);
          });
      stats_update(std::string("rain_flux_sfc"), rain_flux_sfc, stats);
    }

    // Also store the value of surface rain water mixing ratio.

    if (stats.l_sample) {
      Kokkos::parallel_for(
          "advance_microphys::array13", rrm_sfc.extent_int(0), KOKKOS_LAMBDA(int cell) {
            const int index_0 = cell;
            rrm_sfc(index_0)  = rrm_zm(index_0, 0);
          });
      stats_update(std::string("rrm_sfc"), rrm_sfc, stats);
    }

    // iirr > 0
  }

  stats_accumulate_hydromet_api(workspace, gr, ngrdcol, hydromet_dim, hm_metadata, hydromet,
                                rho_ds_zt, stats);

  write_adv_micro_errors(gr, ngrdcol, dt, time_current, hydromet_dim, wm_zt, wp2, exner, rho,
                         rho_zm, rcm, cloud_frac, Kh_zm, Skw_zm, rho_ds_zm, rho_ds_zt,
                         invrs_rho_ds_zt, hydromet_mc, Ncm_mc, Lscale, hydromet_vel_covar_zt_impc,
                         hydromet_vel_covar_zt_expc, clubb_params, nu_vert_res_dep, l_upwind_xm_ma,
                         hydromet, hydromet_vel_zt, hydrometp2, K_hm, Ncm, Nc_in_cloud, rvm_mc,
                         thlm_mc, wphydrometp, wpNcp, err_info);

  return;
}
inline void advance_hydrometeor(
    ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt, int hydromet_dim,
    const hm_metadata_type& hm_metadata, Real2D wm_zt, Real2D exner, Real2D cloud_frac, Real3D K_hm,
    Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real3D hydromet_mc,
    Real3D hydromet_vel_covar_zt_impc, Real3D hydromet_vel_covar_zt_expc,
    const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma, int tridiag_solve_method,
    int fill_holes_type, stats_type& stats, Real3D hydromet, Real3D hydromet_vel_zt,
    Real3D hydrometp2, Real2D rvm_mc, Real2D thlm_mc, err_info_type& err_info, Real3D wphydrometp,
    Real3D hydromet_vel, Real3D hydromet_vel_covar, Real3D hydromet_vel_covar_zt,
    const parameters_microphys& microphys) {
  // Description:
  //   Advance each hydrometeor (precipitating hydrometeor) one model time step.
  // References:
  //   None
  //-----------------------------------------------------------------------
  //---------------------------- Input Variables ----------------------------
  // approximation rather than a centered differencing for turbulent or
  // mean advection terms. It affects rtm, thlm, sclrm, um and vm.
  //---------------------------- Output Variables ----------------------------
  //---------------------------- Local Variables ----------------------------
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope                      = workspace.scope();
  const auto l_hydromet_sed                 = microphys.l_hydromet_sed;
  const auto hm_metadata_hydromet_list_host = to_host_space(workspace, hm_metadata.hydromet_list);
  const bool l_sample                       = stats.l_sample;

  Real3D lhs_ta = workspace_scope.view<Real3D>(
      "advance_hydrometeor::lhs_ta", 3, ngrdcol,
      gr.nzt); // LHS corresponding to contribution from turbulent adv.  [1/s]
  Real3D sed_diff_lhs =
      workspace_scope.view<Real3D>("advance_hydrometeor::sed_diff_lhs", 3, ngrdcol, gr.nzt); // In
  Real3D sed_turb_lhs =
      workspace_scope.view<Real3D>("advance_hydrometeor::sed_turb_lhs", 3, ngrdcol, gr.nzt);
  Real3D lhs_ma = workspace_scope.view<Real3D>(
      "advance_hydrometeor::lhs_ma", 3, ngrdcol,
      gr.nzt); // LHS corresponding to contribution from mean advection  [1/s]
  Real3D lhs = workspace_scope.view<Real3D>("advance_hydrometeor::lhs", 3, ngrdcol,
                                            gr.nzt); // Left hand side array
  Real2D rhs = workspace_scope.view<Real2D>("advance_hydrometeor::rhs", ngrdcol,
                                            gr.nzt); // Right hand side vector
  std::string hydromet_name;
  real_t max_velocity; // Maximum sedimentation velocity         [m/s]
  Real2D hydromet_zm =
      workspace_scope.view<Real2D>("advance_hydrometeor::hydromet_zm", ngrdcol,
                                   gr.nzm); // Mean of hydrometeor interp. to m-levs. [units vary]
  Real3D ratio_hmp2_on_hmm2 =
      workspace_scope.view<Real3D>("advance_hydrometeor::ratio_hmp2_on_hmm2", ngrdcol, gr.nzm,
                                   hydromet_dim); // Value of <hm'^2> / <hm>^2    [-]
  Real2D xpwp       = workspace_scope.view<Real2D>("advance_hydrometeor::xpwp", ngrdcol,
                                                   gr.nzm); // x'w' for arbitrary x
  Real2D K_hm_nu_hm = workspace_scope.view<Real2D>("advance_hydrometeor::K_hm_nu_hm", ngrdcol,
                                                   gr.nzm); // K_hm + nu_vert_res_dep%nu_hm
  constexpr bool l_fill_holes_hm = true;
  int i;   // Loop indices
  int k;   // Loop indices
  int ihm; // Loop indices
  std::string name_bt;
  std::string name_hf;
  std::string name_wvhf;
  std::string name_cl;
  std::string name_mc;
  std::string name_hmp2;
  std::string name_wp;
  const auto gr_nzm                   = gr.nzm;
  const auto hm_metadata_hydromet_tol = hm_metadata.hydromet_tol;
  const auto nu_vert_res_dep_nu_hm    = nu_vert_res_dep.nu_hm;

  //---------------------------- Begin Code ----------------------------

  // Loop over each type of precipitating hydrometeor and advance one model
  // time step.
  for (ihm = 0; ihm < hydromet_dim; ++ihm) {

    // Set up hydrometeor stats names and velocity caps.
    setup_stats_names(ihm, hydromet_dim, hm_metadata_hydromet_list_host, max_velocity, &name_bt,
                      &name_hf, &name_wvhf, &name_cl, &name_mc);

    if (stats.l_sample) {
      if ((static_cast<int>(microphys_detail::trim(name_mc).size()) > 0)) {
        // Explicit microphysics contribution for this hydrometeor species.
        // Update explicit contributions to the hydrometeor species.
        stats_update(name_mc, Kokkos::subview(hydromet_mc, Kokkos::ALL, Kokkos::ALL, ihm), stats);
      }
      if ((static_cast<int>(microphys_detail::trim(name_bt).size()) > 0)) {
        // Save prior value of the hydrometeors for determining total time
        // tendency.
        Real2D array_actual_1 = workspace_scope.view<Real2D>(
            "advance_hydrometeor::array_actual_1", hydromet.extent_int(0), hydromet.extent_int(1));
        Kokkos::parallel_for(
            "advance_hydrometeor::array2", ((hydromet.extent_int(0) * hydromet.extent_int(1))),
            KOKKOS_LAMBDA(int cell) {
              const int index_0                = ((cell % hydromet.extent_int(0)));
              const int index_1                = ((cell / hydromet.extent_int(0)));
              array_actual_1(index_0, index_1) = (hydromet(index_0, index_1, ihm) / dt);
            });
        stats_begin_budget(name_bt, array_actual_1, stats);
      }
    }

    // Set realistic limits on sedimentation velocities, following the
    // numbers in the Morrison microphysics.
    for (k = 0; k < gr.nzt; ++k) {
      for (i = 0; i < ngrdcol; ++i) {
        if (clubb_at_least_debug_level_api(1)) {
          // Print a warning if the velocity has a large magnitude or the
          // velocity is in the wrong direction.
          if (((to_host_space(workspace, hydromet_vel_zt)(i, k, ihm) < max_velocity) ||
               (to_host_space(workspace, hydromet_vel_zt)(i, k, ihm) > zero_threshold))) {

            std::cerr << (microphys_detail::trim(
                              std::string(hm_metadata_hydromet_list_host.data() +
                                              ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                          hm_metadata_hydromet_list_host.extent(0))) +
                          std::string(" velocity at i, k = "))
                      << " " << (((i + 1))) << " " << (((k + 1))) << " " << std::string(" = ")
                      << " " << to_host_space(workspace, hydromet_vel_zt)(i, k, ihm) << " "
                      << std::string("m/s") << "\n";
          }
        }
        Kokkos::parallel_for(
            "advance_hydrometeor::element", 1, KOKKOS_LAMBDA(int) {
              hydromet_vel_zt(i, k, ihm) = ((((hydromet_vel_zt(i, k, ihm)) > (max_velocity)
                                                  ? (hydromet_vel_zt(i, k, ihm))
                                                  : (max_velocity))) < (zero_threshold)
                                                ? (((hydromet_vel_zt(i, k, ihm)) > (max_velocity)
                                                        ? (hydromet_vel_zt(i, k, ihm))
                                                        : (max_velocity)))
                                                : (zero_threshold));
            });
      }
      // Native loop: k = 1, gr%nzt, 1
    }

    // Interpolate velocity to the momentum grid for a centered difference
    // approximation of the sedimenation term.
    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr,
              Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, ihm),
              Kokkos::subview(hydromet_vel, Kokkos::ALL, Kokkos::ALL, ihm));
    // Upper boundary condition
    Kokkos::parallel_for(
        "advance_hydrometeor::array3", hydromet_vel.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0                          = cell;
          hydromet_vel(index_0, ((gr_nzm - 1)), ihm) = zero;
        });

    // Calculate the value of <hm'^2> / <hm>^2.  This will be used to update
    // <hm'^2> after <hm> has been advanced one model timestep.  This method
    // is being used because CLUBB does not currently have a predictive
    // equation for <hm'^2> (hydrometp2).
    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm),
              hydromet_zm);

    Kokkos::parallel_for(
        "advance_hydrometeor::loop4", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzm; ++k) {
            if ((hydromet_zm(i, k) > hm_metadata_hydromet_tol(ihm))) {

              // Calculate the ratio of the overall variance of the hydrometeor
              // to the overall mean of the hydrometeor squared.
              ratio_hmp2_on_hmm2(i, k, ihm) =
                  (hydrometp2(i, k, ihm) / (hydromet_zm(i, k) * hydromet_zm(i, k)));

            } else {
              // hydromet_zm(i,k) <= hydromet_tol(ihm)

              // The overall mean of the hydrometeor is 0 or is treated as 0.  The
              // overall variance of the hydrometeor is also 0 in this scenario.
              // The ratio is undefined, but will be assigned a value of 0.
              ratio_hmp2_on_hmm2(i, k, ihm) = zero;

              // hydromet_zm(i,k) > hydromet_tol(ihm)
            }
            // Native loop: k = 1, gr%nzm, 1
          }
        });

    // Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
    // a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
    // A Crank-Nicholson time-stepping scheme is used for this variable.
    // This is the portion of the calculation using < h_m > from timestep t.
    Kokkos::parallel_for(
        "advance_hydrometeor::loop5", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzm; ++k) {
            K_hm_nu_hm(i, k) = (K_hm(i, k, ihm) + nu_vert_res_dep_nu_hm(i));
          }
        });

    calc_xpwp(gr.nzm, gr.nzt, ngrdcol, gr, K_hm_nu_hm,
              Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm), xpwp);

    Kokkos::parallel_for(
        "advance_hydrometeor::array6",
        ((wphydrometp.extent_int(0) * gr.nzm - 2 * wphydrometp.extent_int(0))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % wphydrometp.extent_int(0)));
          const int index_1 = ((cell / wphydrometp.extent_int(0)));
          wphydrometp(index_0, ((index_1 + 1)), ihm) =
              ((-one_half) * xpwp(index_0, ((index_1 + 1))));
        });

    // A zero-flux boundary condition is used for hydrometeors.
    Kokkos::parallel_for(
        "advance_hydrometeor::array7", wphydrometp.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0            = cell;
          wphydrometp(index_0, 0, ihm) = zero;
        });
    Kokkos::parallel_for(
        "advance_hydrometeor::array8", wphydrometp.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0                         = cell;
          wphydrometp(index_0, ((gr_nzm - 1)), ihm) = zero;
        });

    // Add implicit terms to the LHS matrix
    microphys_lhs(
        workspace, gr, ngrdcol,
        microphys_detail::trim(std::string(hm_metadata_hydromet_list_host.data() +
                                               ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                           hm_metadata_hydromet_list_host.extent(0))),
        to_host_space(workspace, l_hydromet_sed)(ihm), dt,
        Kokkos::subview(K_hm, Kokkos::ALL, Kokkos::ALL, ihm), nu_vert_res_dep.nu_hm, wm_zt,
        Kokkos::subview(hydromet_vel, Kokkos::ALL, Kokkos::ALL, ihm),
        Kokkos::subview(hydromet_vel_zt, Kokkos::ALL, Kokkos::ALL, ihm),
        Kokkos::subview(hydromet_vel_covar_zt_impc, Kokkos::ALL, Kokkos::ALL, ihm), rho_ds_zm,
        rho_ds_zt, invrs_rho_ds_zt, l_upwind_xm_ma, stats, lhs_ta, lhs_ma, sed_turb_lhs,
        sed_diff_lhs, lhs, microphys);

    // Set up explicit term in the RHS vector
    microphys_rhs(
        workspace, gr, ngrdcol,
        microphys_detail::trim(std::string(hm_metadata_hydromet_list_host.data() +
                                               ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                           hm_metadata_hydromet_list_host.extent(0))),
        dt, to_host_space(workspace, l_hydromet_sed)(ihm),
        Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm),
        Kokkos::subview(hydromet_mc, Kokkos::ALL, Kokkos::ALL, ihm),
        Kokkos::subview(K_hm, Kokkos::ALL, Kokkos::ALL, ihm), nu_vert_res_dep.nu_hm, cloud_frac,
        Kokkos::subview(hydromet_vel_covar_zt_expc, Kokkos::ALL, Kokkos::ALL, ihm), rho_ds_zm,
        rho_ds_zt, invrs_rho_ds_zt, stats, rhs, microphys);

    //!!!! Advance hydrometeor one time step.
    microphys_solve(
        workspace, gr, ngrdcol,
        microphys_detail::trim(std::string(hm_metadata_hydromet_list_host.data() +
                                               ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                           hm_metadata_hydromet_list_host.extent(0))),
        to_host_space(workspace, l_hydromet_sed)(ihm), lhs_ta, lhs_ma, sed_turb_lhs, sed_diff_lhs,
        cloud_frac, tridiag_solve_method, stats, lhs, rhs,
        Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm), err_info, microphys);

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << "\n";
        std::cerr << (std::string("Error in hydrometeor field ") +
                      microphys_detail::trim(
                          std::string(hm_metadata_hydromet_list_host.data() +
                                          ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                      hm_metadata_hydromet_list_host.extent(0))))
                  << "\n";

        std::cerr << (microphys_detail::trim(
                          std::string(hm_metadata_hydromet_list_host.data() +
                                          ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                      hm_metadata_hydromet_list_host.extent(0))) +
                      std::string(" = "))
                  << " "
                  << microphys_detail::dump_values(
                         Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm))
                  << "\n";

        // err_info%err_code == clubb_fatal_error
      }
    }

    // Native loop: ihm = 1, hydromet_dim, 1
  }

  // Now that all precipitating hydrometeors have been advanced, fill holes in
  // hydromet profiles.
  fill_holes_driver_api(workspace, gr, ngrdcol, gr.nzt, dt, hydromet_dim, hm_metadata,
                        l_fill_holes_hm, rho_ds_zt, exner, fill_holes_type, stats, thlm_mc, rvm_mc,
                        hydromet);

  // Loop over each type of precipitating hydrometeor and calculate hydrometeor
  // covariances (<w'hm'> and <V_hm'hm'>) and other quantities requiring the
  // value of hydromet (<hm>) from the (t+1) timestep.
  for (ihm = 0; ihm < hydromet_dim; ++ihm) {

    // Set up hydrometeor stats names and velocity caps.
    setup_stats_names(ihm, hydromet_dim, hm_metadata_hydromet_list_host, max_velocity, &name_bt,
                      &name_hf, &name_wvhf, &name_cl, &name_mc);

    // Print warning message if any hydrometeor species has a value < 0.
    if (microphys_detail::any_less(Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm),
                                   zero_threshold)) {

      hydromet_name = std::string(hm_metadata_hydromet_list_host.data() +
                                      ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                  hm_metadata_hydromet_list_host.extent(0));

      if (clubb_at_least_debug_level_api(1)) {
        for (k = 0; k < gr.nzt; ++k) {
          for (i = 0; i < ngrdcol; ++i) {
            if ((to_host_space(workspace, hydromet)(i, k, ihm) < zero_threshold)) {
              std::cerr << (microphys_detail::trim(hydromet_name) + std::string(" < ")) << " "
                        << zero_threshold << " " << std::string(" in advance_microphys at i, k= ")
                        << " " << (((i + 1))) << " " << (((k + 1))) << "\n";
            }
          }
        }
      }

      // hydromet(:,:,ihm) < 0
    }

    // Lower boundary condition
    // Hydrometeors that are below the model lower boundary level have
    // sedimented out of the model domain, and is not conserved.
    Kokkos::parallel_for(
        "advance_hydrometeor::loop9", ngrdcol, KOKKOS_LAMBDA(int column) {
          const int i = column;
          if ((hydromet(i, 0, ihm) < hm_metadata_hydromet_tol(ihm))) {
            hydromet(i, 0, ihm) = zero_threshold;
          }
        });

    // Calculate the value of <hm'^2> (hydrometp2).  This is based on the
    // ratio of <hm'^2> / <hm>^2 (ratio_hmp2_on_hmm2) that was saved before
    // hydrometeors were updated.  This method is being used because CLUBB
    // does not currently have a predictive equation for <hm'^2>.
    auto interpolated_11 = zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr,
                                     Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm));
    Kokkos::parallel_for(
        "advance_hydrometeor::array11", ((hydromet_zm.extent_int(0) * hydromet_zm.extent_int(1))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0             = ((cell % hydromet_zm.extent_int(0)));
          const int index_1             = ((cell / hydromet_zm.extent_int(0)));
          hydromet_zm(index_0, index_1) = ((interpolated_11(index_0, index_1)) > (real_t(0.0))
                                               ? (interpolated_11(index_0, index_1))
                                               : (real_t(0.0)));
        });

    Kokkos::parallel_for(
        "advance_hydrometeor::loop12", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzm; ++k) {
            hydrometp2(i, k, ihm) =
                (ratio_hmp2_on_hmm2(i, k, ihm) * (hydromet_zm(i, k) * hydromet_zm(i, k)));
            // Native loop: k = 1, gr%nz, 1
          }
        });

    // Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
    // a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
    // A Crank-Nicholson time-stepping scheme is used for this variable.
    // This is the portion of the calculation using < h_m > from timestep t+1.
    Kokkos::parallel_for(
        "advance_hydrometeor::loop13", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzm; ++k) {
            K_hm_nu_hm(i, k) = (K_hm(i, k, ihm) + nu_vert_res_dep_nu_hm(i));
          }
        });

    calc_xpwp(gr.nzm, gr.nzt, ngrdcol, gr, K_hm_nu_hm,
              Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, ihm), xpwp);

    Kokkos::parallel_for(
        "advance_hydrometeor::array14",
        ((wphydrometp.extent_int(0) * gr.nzm - 2 * wphydrometp.extent_int(0))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % wphydrometp.extent_int(0)));
          const int index_1 = ((cell / wphydrometp.extent_int(0)));
          wphydrometp(index_0, ((index_1 + 1)), ihm) =
              ((-one_half) * xpwp(index_0, ((index_1 + 1))));
        });

    // A zero-flux boundary condition is used for hydrometeors.
    Kokkos::parallel_for(
        "advance_hydrometeor::array15", wphydrometp.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0            = cell;
          wphydrometp(index_0, 0, ihm) = zero;
        });
    Kokkos::parallel_for(
        "advance_hydrometeor::array16", wphydrometp.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0                         = cell;
          wphydrometp(index_0, ((gr_nzm - 1)), ihm) = zero;
        });

    //!! Calculate the covariance of hydrometeor sedimentation velocity and
    //!! the hydrometeor, which is solved semi-implicitly on thermodynamic
    //!! levels.
    Kokkos::parallel_for(
        "advance_hydrometeor::array17",
        ((hydromet_vel_covar_zt.extent_int(0) * hydromet_vel_covar_zt.extent_int(1))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % hydromet_vel_covar_zt.extent_int(0)));
          const int index_1 = ((cell / hydromet_vel_covar_zt.extent_int(0)));
          hydromet_vel_covar_zt(index_0, index_1, ihm) =
              ((hydromet_vel_covar_zt_impc(index_0, index_1, ihm) *
                hydromet(index_0, index_1, ihm)) +
               hydromet_vel_covar_zt_expc(index_0, index_1, ihm));
        });

    //!! Calculate the covariance of hydrometeor sedimentation velocity and
    //!! the hydrometeor, < V_hm'h_m' >, by interpolating the thermodynamic
    //!! level results to momentum levels.
    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr,
              Kokkos::subview(hydromet_vel_covar_zt, Kokkos::ALL, Kokkos::ALL, ihm),
              Kokkos::subview(hydromet_vel_covar, Kokkos::ALL, Kokkos::ALL, ihm));

    // Boundary conditions for < V_hm'hm' >.
    Kokkos::parallel_for(
        "advance_hydrometeor::array18", hydromet_vel_covar.extent_int(0), KOKKOS_LAMBDA(int cell) {
          const int index_0                                = cell;
          hydromet_vel_covar(index_0, ((gr_nzm - 1)), ihm) = zero;
        });

    // Statistics for all covariances involving hydrometeors:  < w'h_m' >,
    // <V_rr'r_r'>, and <V_Nr'N_r'>.
    // Covariance of vertical velocity and the hydrometeor.
    // Covariance of sedimentation velocity of r_r and r_r.
    // Covariance of sedimentation velocity of N_r and N_r.

    if (stats.l_sample) {
      if ((static_cast<int>(microphys_detail::trim(
                                std::string(hm_metadata_hydromet_list_host.data() +
                                                ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                            hm_metadata_hydromet_list_host.extent(0)))
                                .size()) > 0)) {
        name_hmp2 = (microphys_detail::trim(
                         std::string(hm_metadata_hydromet_list_host.data() +
                                         ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                     hm_metadata_hydromet_list_host.extent(0))
                             .substr(0, 2)) +
                     std::string("p2"));
        name_wp   = ((std::string("wp") +
                      microphys_detail::trim(
                          std::string(hm_metadata_hydromet_list_host.data() +
                                          ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                                      hm_metadata_hydromet_list_host.extent(0))
                              .substr(0, 2))) +
                     std::string("p"));
        stats_update(name_hmp2, Kokkos::subview(hydrometp2, Kokkos::ALL, Kokkos::ALL, ihm), stats);
        // This section has not been GPUized yet, but this variable has a GPU copy.
        // Update that copy with the CPU result before stats reads it back.
        stats_update(name_wp, Kokkos::subview(wphydrometp, Kokkos::ALL, Kokkos::ALL, ihm), stats);
      }
      if ((microphys_detail::trim(
               std::string(hm_metadata_hydromet_list_host.data() +
                               ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                           hm_metadata_hydromet_list_host.extent(0))) == std::string("rrm"))) {
        stats_update(
            std::string("Vrrprrp"),
            Kokkos::subview(hydromet_vel_covar, Kokkos::ALL, Kokkos::ALL, hm_metadata.iirr), stats);
      } else {
        if ((microphys_detail::trim(
                 std::string(hm_metadata_hydromet_list_host.data() +
                                 ((hm_metadata_hydromet_list_host.extent(0) * ihm)),
                             hm_metadata_hydromet_list_host.extent(0))) == std::string("Nrm"))) {
          stats_update(
              std::string("VNrpNrp"),
              Kokkos::subview(hydromet_vel_covar, Kokkos::ALL, Kokkos::ALL, hm_metadata.iiNr),
              stats);
        }
      }
    }

    if (stats.l_sample) {
      if ((static_cast<int>(microphys_detail::trim(name_bt).size()) > 0)) {
        // Total time tendency
        Real2D array_actual_19 = workspace_scope.view<Real2D>(
            "advance_hydrometeor::array_actual_19", hydromet.extent_int(0), hydromet.extent_int(1));
        Kokkos::parallel_for(
            "advance_hydrometeor::array20", ((hydromet.extent_int(0) * hydromet.extent_int(1))),
            KOKKOS_LAMBDA(int cell) {
              const int index_0                 = ((cell % hydromet.extent_int(0)));
              const int index_1                 = ((cell / hydromet.extent_int(0)));
              array_actual_19(index_0, index_1) = (hydromet(index_0, index_1, ihm) / dt);
            });
        stats_finalize_budget(name_bt, array_actual_19, stats);
      }
    }

    // Native loop: ihm = 1, hydromet_dim, 1
  }

  return;
}
inline void advance_Ncm(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                        Real2D wm_zt, Real2D cloud_frac, Real2D K_Nc, Real2D rcm, Real2D rho_ds_zm,
                        Real2D rho_ds_zt, Real2D invrs_rho_ds_zt, Real2D Ncm_mc,
                        const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma,
                        int tridiag_solve_method, stats_type& stats, Real2D Ncm, Real2D Nc_in_cloud,
                        err_info_type& err_info, Real2D wpNcp,
                        const parameters_microphys& microphys) {
  // Description:
  // Advance cloud droplet concentration (Ncm) one model time step.
  // References:
  //-----------------------------------------------------------------------
  //------------------------------- Input Variables -------------------------------
  // approximation rather than a centered differencing for turbulent or
  // mean advection terms. It affects rtm, thlm, sclrm, um and vm.
  //------------------------------- Output Variables -------------------------------
  //------------------------------- Locak Variables -------------------------------
  // Sedimentation velocity of cloud droplet concentration is not considered in
  // the solution of N_c.
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope          = workspace.scope();
  const auto l_in_cloud_Nc_diff = microphys.l_in_cloud_Nc_diff;
  const bool l_sample           = stats.l_sample;

  Real2D Ncm_vel_covar_zt_impc =
      workspace_scope.view<Real2D>("advance_Ncm::Ncm_vel_covar_zt_impc", ngrdcol, gr.nzt); // In
  Real2D Ncm_vel_covar_zt_expc =
      workspace_scope.view<Real2D>("advance_Ncm::Ncm_vel_covar_zt_expc", ngrdcol,
                                   gr.nzt); // Exp. comp. <V_Nc'N_c'> t-levs [(num/kg)(m/s)]
  Real2D Ncm_vel = workspace_scope.view<Real2D>(
      "advance_Ncm::Ncm_vel", ngrdcol,
      gr.nzm); // Mean N_c sedimentation velocity, <V_Nc> (m-levs.) [m/s]
  Real2D Ncm_vel_zt =
      workspace_scope.view<Real2D>("advance_Ncm::Ncm_vel_zt", ngrdcol, gr.nzt); // In
  Real3D lhs =
      workspace_scope.view<Real3D>("advance_Ncm::lhs", 3, ngrdcol, gr.nzt); // Left hand side array
  Real2D rhs =
      workspace_scope.view<Real2D>("advance_Ncm::rhs", ngrdcol, gr.nzt); // Right hand side vector
  Real2D Ncm_mvr_min = workspace_scope.view<Real2D>(
      "advance_Ncm::Ncm_mvr_min", ngrdcol,
      gr.nzt); // Min. allowable Ncm due to max. cloud droplet mvr [num/kg]
  Real2D Ncm_min = workspace_scope.view<Real2D>(
      "advance_Ncm::Ncm_min", ngrdcol,
      gr.nzt); // Minimum allowable mean cloud droplet conc.       [num/kg]
  Real2D Ncic_min = workspace_scope.view<Real2D>(
      "advance_Ncm::Ncic_min", ngrdcol,
      gr.nzt); // Min. allowable in-cloud mean cloud droplet conc. [num/kg]
  Real2D xpwp =
      workspace_scope.view<Real2D>("advance_Ncm::xpwp", ngrdcol, gr.nzm); // x'w' for arbitrary x
  Real2D K_Nc_nu_hm = workspace_scope.view<Real2D>("advance_Ncm::K_Nc_nu_hm", ngrdcol,
                                                   gr.nzm); // K_Nc + nu_vert_res_dep%nu_hm
  Real3D lhs_ta     = workspace_scope.view<Real3D>(
      "advance_Ncm::lhs_ta", 3, ngrdcol,
      gr.nzt); // LHS corresponding to contribution from turbulent adv.  [1/s]
  Real3D sed_diff_lhs =
      workspace_scope.view<Real3D>("advance_Ncm::sed_diff_lhs", 3, ngrdcol, gr.nzt); // In
  Real3D sed_turb_lhs =
      workspace_scope.view<Real3D>("advance_Ncm::sed_turb_lhs", 3, ngrdcol, gr.nzt);
  Real3D lhs_ma = workspace_scope.view<Real3D>(
      "advance_Ncm::lhs_ma", 3, ngrdcol,
      gr.nzt); // LHS corresponding to contribution from mean advection  [1/s]
  constexpr bool l_Ncm_sed = false;
  int i; // Loop indices
  int k; // Loop indices
  const auto gr_nzm                = gr.nzm;
  const auto gr_nzt                = gr.nzt;
  const auto nu_vert_res_dep_nu_hm = nu_vert_res_dep.nu_hm;

  //------------------------------- Begin Code -------------------------------

  // The mean sedimentation velocity of cloud droplet concentration, < V_Nc >,
  // and the covariance of N_c sedimentation velocity with N_c, < V_Nc'Nc' >,
  // are not used in the advancement of < N_c >.  Sedimentation velocity of N_c
  // is not considered in the solution.  However, these variables need to be
  // initialized to 0 and passed into microphys_lhs and microphys_rhs.
  Kokkos::parallel_for(
      "advance_ncm::array1",
      ((Ncm_vel_covar_zt_impc.extent_int(0) * Ncm_vel_covar_zt_impc.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0                       = ((cell % Ncm_vel_covar_zt_impc.extent_int(0)));
        const int index_1                       = ((cell / Ncm_vel_covar_zt_impc.extent_int(0)));
        Ncm_vel_covar_zt_impc(index_0, index_1) = zero;
      });
  Kokkos::parallel_for(
      "advance_ncm::array2",
      ((Ncm_vel_covar_zt_expc.extent_int(0) * Ncm_vel_covar_zt_expc.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0                       = ((cell % Ncm_vel_covar_zt_expc.extent_int(0)));
        const int index_1                       = ((cell / Ncm_vel_covar_zt_expc.extent_int(0)));
        Ncm_vel_covar_zt_expc(index_0, index_1) = zero;
      });
  Kokkos::parallel_for(
      "advance_ncm::array3", ((Ncm_vel.extent_int(0) * Ncm_vel.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0         = ((cell % Ncm_vel.extent_int(0)));
        const int index_1         = ((cell / Ncm_vel.extent_int(0)));
        Ncm_vel(index_0, index_1) = zero;
      });
  Kokkos::parallel_for(
      "advance_ncm::array4", ((Ncm_vel_zt.extent_int(0) * Ncm_vel_zt.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0            = ((cell % Ncm_vel_zt.extent_int(0)));
        const int index_1            = ((cell / Ncm_vel_zt.extent_int(0)));
        Ncm_vel_zt(index_0, index_1) = zero;
      });

  if (stats.l_sample) {
    // Update explicit contributions to cloud droplet concentration.
    stats_update(std::string("Ncm_mc"), Ncm_mc, stats);

    if (l_in_cloud_Nc_diff) {
      // Save prior value of Ncm for determining total time tendency.
      Real2D array_actual_5 = workspace_scope.view<Real2D>(
          "advance_ncm::array_actual_5", Nc_in_cloud.extent_int(0), Nc_in_cloud.extent_int(1));
      Kokkos::parallel_for(
          "advance_ncm::array6", ((Nc_in_cloud.extent_int(0) * Nc_in_cloud.extent_int(1))),
          KOKKOS_LAMBDA(int cell) {
            const int index_0 = ((cell % Nc_in_cloud.extent_int(0)));
            const int index_1 = ((cell / Nc_in_cloud.extent_int(0)));
            array_actual_5(index_0, index_1) =
                ((Nc_in_cloud(index_0, index_1) * ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                                       ? (cloud_frac(index_0, index_1))
                                                       : (cloud_frac_min))) /
                 dt);
          });
      stats_begin_budget(std::string("Ncm_bt"), array_actual_5, stats);
    } else {
      // Save prior value of Ncm for determining total time tendency.
      Real2D array_actual_7 = workspace_scope.view<Real2D>("advance_ncm::array_actual_7",
                                                           Ncm.extent_int(0), Ncm.extent_int(1));
      Kokkos::parallel_for(
          "advance_ncm::array8", ((Ncm.extent_int(0) * Ncm.extent_int(1))),
          KOKKOS_LAMBDA(int cell) {
            const int index_0                = ((cell % Ncm.extent_int(0)));
            const int index_1                = ((cell / Ncm.extent_int(0)));
            array_actual_7(index_0, index_1) = (Ncm(index_0, index_1) / dt);
          });
      stats_begin_budget(std::string("Ncm_bt"), array_actual_7, stats);
    }
  }

  // Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
  // a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
  // A Crank-Nicholson time-stepping scheme is used for this variable.
  // This is the portion of the calculation using < N_c > from timestep t.
  Kokkos::parallel_for(
      "advance_ncm::loop9", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Loop indices
        const int i = column;
        for (k = 0; k < gr_nzm; ++k) {
          K_Nc_nu_hm(i, k) = (K_Nc(i, k) + nu_vert_res_dep_nu_hm(i));
        }
      });

  calc_xpwp(gr.nzm, gr.nzt, ngrdcol, gr, K_Nc_nu_hm, Ncm, xpwp);

  Kokkos::parallel_for(
      "advance_ncm::array10", ((wpNcp.extent_int(0) * gr.nzm - 2 * wpNcp.extent_int(0))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0               = ((cell % wpNcp.extent_int(0)));
        const int index_1               = ((cell / wpNcp.extent_int(0)));
        wpNcp(index_0, ((index_1 + 1))) = ((-one_half) * xpwp(index_0, ((index_1 + 1))));
      });

  // A zero-flux boundary condition is used for N_c.
  Kokkos::parallel_for(
      "advance_ncm::array11", wpNcp.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0 = cell;
        wpNcp(index_0, 0) = zero;
      });
  Kokkos::parallel_for(
      "advance_ncm::array12", wpNcp.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0              = cell;
        wpNcp(index_0, ((gr_nzm - 1))) = zero;
      });

  // Covariance of vertical velocity and N_c.

  // Add implicit terms to the LHS array
  microphys_lhs(workspace, gr, ngrdcol, std::string("Ncm"), l_Ncm_sed, dt, K_Nc,
                nu_vert_res_dep.nu_hm, wm_zt, Ncm_vel, Ncm_vel_zt, Ncm_vel_covar_zt_impc, rho_ds_zm,
                rho_ds_zt, invrs_rho_ds_zt, l_upwind_xm_ma, stats, lhs_ta, lhs_ma, sed_turb_lhs,
                sed_diff_lhs, lhs, microphys);

  // Set up explicit term in the RHS vector
  if (l_in_cloud_Nc_diff) {

    Real2D array_actual_13 = workspace_scope.view<Real2D>(
        "advance_ncm::array_actual_13", Ncm_mc.extent_int(0), Ncm_mc.extent_int(1));
    Kokkos::parallel_for(
        "advance_ncm::array14", ((Ncm_mc.extent_int(0) * Ncm_mc.extent_int(1))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % Ncm_mc.extent_int(0)));
          const int index_1 = ((cell / Ncm_mc.extent_int(0)));
          array_actual_13(index_0, index_1) =
              (Ncm_mc(index_0, index_1) / ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                               ? (cloud_frac(index_0, index_1))
                                               : (cloud_frac_min)));
        });
    microphys_rhs(workspace, gr, ngrdcol, std::string("Ncm"), dt, l_Ncm_sed, Nc_in_cloud,
                  array_actual_13, K_Nc, nu_vert_res_dep.nu_hm, cloud_frac, Ncm_vel_covar_zt_expc,
                  rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, stats, rhs, microphys);

  } else {

    microphys_rhs(workspace, gr, ngrdcol, std::string("Ncm"), dt, l_Ncm_sed, Ncm, Ncm_mc, K_Nc,
                  nu_vert_res_dep.nu_hm, cloud_frac, Ncm_vel_covar_zt_expc, rho_ds_zm, rho_ds_zt,
                  invrs_rho_ds_zt, stats, rhs, microphys);
  }

  //!!!! Advance Ncm one time step.
  if (l_in_cloud_Nc_diff) {

    microphys_solve(workspace, gr, ngrdcol, std::string("Ncm"), l_Ncm_sed, lhs_ta, lhs_ma,
                    sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method, stats, lhs, rhs,
                    Nc_in_cloud, err_info, microphys);

    Kokkos::parallel_for(
        "advance_ncm::array15", ((Ncm.extent_int(0) * Ncm.extent_int(1))), KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % Ncm.extent_int(0)));
          const int index_1 = ((cell / Ncm.extent_int(0)));
          Ncm(index_0, index_1) =
              (Nc_in_cloud(index_0, index_1) * ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                                    ? (cloud_frac(index_0, index_1))
                                                    : (cloud_frac_min)));
        });

  } else {

    microphys_solve(workspace, gr, ngrdcol, std::string("Ncm"), l_Ncm_sed, lhs_ta, lhs_ma,
                    sed_turb_lhs, sed_diff_lhs, cloud_frac, tridiag_solve_method, stats, lhs, rhs,
                    Ncm, err_info, microphys);

    Kokkos::parallel_for(
        "advance_ncm::array16", ((Nc_in_cloud.extent_int(0) * Nc_in_cloud.extent_int(1))),
        KOKKOS_LAMBDA(int cell) {
          const int index_0 = ((cell % Nc_in_cloud.extent_int(0)));
          const int index_1 = ((cell / Nc_in_cloud.extent_int(0)));
          Nc_in_cloud(index_0, index_1) =
              (Ncm(index_0, index_1) / ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                            ? (cloud_frac(index_0, index_1))
                                            : (cloud_frac_min)));
        });
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << std::string("calling microphys_solve in advance_Ncm") << "\n";
      return;
    }
  }

  // Clipping for mean cloud droplet concentration, <Nc>.

  // When mean cloud water mixing ratio, <rc>, is found at a grid level, mean
  // cloud droplet concentration must be at least a minimum value so that
  // average cloud droplet mean volume radius stays within an upper bound.

  // The minimum mean rain drop concentration is given by:
  //
  // <Nc> = <rc> / ( (4/3) * pi * rho_lw * mvr_cloud_max^3 ).

  Kokkos::parallel_for(
      "advance_ncm::array17", ((Ncm_mvr_min.extent_int(0) * Ncm_mvr_min.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % Ncm_mvr_min.extent_int(0)));
        const int index_1 = ((cell / Ncm_mvr_min.extent_int(0)));
        Ncm_mvr_min(index_0, index_1) =
            ((one /
              (((four_thirds * pi) * rho_lw) * ((mvr_cloud_max * mvr_cloud_max) * mvr_cloud_max))) *
             rcm(index_0, index_1));
      });

  // Calculate the minimum threshold value for mean cloud droplet
  // concentration.  This is the greater of the minimum allowable Ncm based on
  // the maximum cloud droplet mean volume radius and Nc_in_cloud_min
  // multiplied by max( cloud_frac, cloud_frac_min ).
  Kokkos::parallel_for(
      "advance_ncm::array18", ((Ncm_min.extent_int(0) * Ncm_min.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % Ncm_min.extent_int(0)));
        const int index_1 = ((cell / Ncm_min.extent_int(0)));
        Ncm_min(index_0, index_1) =
            (((Nc_in_cloud_min * ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                      ? (cloud_frac(index_0, index_1))
                                      : (cloud_frac_min)))) > (Ncm_mvr_min(index_0, index_1))
                 ? ((Nc_in_cloud_min * ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                            ? (cloud_frac(index_0, index_1))
                                            : (cloud_frac_min))))
                 : (Ncm_mvr_min(index_0, index_1)));
      });

  // Print warning message if Ncm (at any level) has a value < Ncm_min.
  if (microphys_detail::any_less(Ncm, Ncm_min)) {

    if (clubb_at_least_debug_level_api(1)) {
      for (k = 0; k < gr.nzt; ++k) {
        for (i = 0; i < ngrdcol; ++i) {
          if ((to_host_space(workspace, Ncm)(i, k) < to_host_space(workspace, Ncm_min)(i, k))) {
            std::cerr << std::string("Ncm < ") << " " << to_host_space(workspace, Ncm_min)(i, k)
                      << " " << std::string(" in advance_microphys at i, k = ") << " "
                      << (((i + 1))) << " " << (((k + 1))) << "\n";
            // Ncm(i,k) < Ncm_min(i,k)
          }
        }
        // Native loop: k = 1, gr%nzt, 1
      }
      // clubb_at_least_debug_level_api( 1 )
    }
    // Ncm < Ncm_min
  }
  // Store the previous value of Ncm for the effect of clipping.

  if (stats.l_sample) {
    Real2D array_actual_19 = workspace_scope.view<Real2D>("advance_ncm::array_actual_19",
                                                          Ncm.extent_int(0), Ncm.extent_int(1));
    Kokkos::parallel_for(
        "advance_ncm::array20", ((Ncm.extent_int(0) * Ncm.extent_int(1))), KOKKOS_LAMBDA(int cell) {
          const int index_0                 = ((cell % Ncm.extent_int(0)));
          const int index_1                 = ((cell / Ncm.extent_int(0)));
          array_actual_19(index_0, index_1) = (Ncm(index_0, index_1) / dt);
        });
    stats_begin_budget(std::string("Ncm_cl"), array_actual_19, stats);
  }

  if (microphys_detail::any_less(Ncm, Ncm_min)) {

    // Clip any values of Ncm below the minimum threshold value of Ncm_min to
    // the minimum threshold value.
    Kokkos::parallel_for(
        "advance_ncm::loop21", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i = column;
          for (k = 0; k < gr_nzt; ++k) {
            if ((Ncm(i, k) < Ncm_min(i, k))) {
              Ncm(i, k) = Ncm_min(i, k);
            }
          }
        });

    // Ncm < Ncm_min
  }
  // Enter the new value of Ncm for the effect of clipping.

  if (stats.l_sample) {
    Real2D array_actual_22 = workspace_scope.view<Real2D>("advance_ncm::array_actual_22",
                                                          Ncm.extent_int(0), Ncm.extent_int(1));
    Kokkos::parallel_for(
        "advance_ncm::array23", ((Ncm.extent_int(0) * Ncm.extent_int(1))), KOKKOS_LAMBDA(int cell) {
          const int index_0                 = ((cell % Ncm.extent_int(0)));
          const int index_1                 = ((cell / Ncm.extent_int(0)));
          array_actual_22(index_0, index_1) = (Ncm(index_0, index_1) / dt);
        });
    stats_finalize_budget(std::string("Ncm_cl"), array_actual_22, stats);
  }

  // Calculate the minimum threshold value for in-cloud mean cloud droplet
  // concentration.  This is the greater of the minimum allowable Nc_in_cloud
  // based on the maximum cloud droplet mean volume radius (Ncm_mvr_min)
  // divided by max( cloud_frac, cloud_frac_min ) and Nc_in_cloud_min.
  Kokkos::parallel_for(
      "advance_ncm::array24", ((Ncic_min.extent_int(0) * Ncic_min.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % Ncic_min.extent_int(0)));
        const int index_1 = ((cell / Ncic_min.extent_int(0)));
        Ncic_min(index_0, index_1) =
            ((Nc_in_cloud_min) > ((Ncm_mvr_min(index_0, index_1) /
                                   ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                                        ? (cloud_frac(index_0, index_1))
                                        : (cloud_frac_min))))
                 ? (Nc_in_cloud_min)
                 : ((Ncm_mvr_min(index_0, index_1) /
                     ((cloud_frac(index_0, index_1)) > (cloud_frac_min)
                          ? (cloud_frac(index_0, index_1))
                          : (cloud_frac_min)))));
      });

  // Clip any negative values of Nc_in_cloud to the minimum threshold value of
  // Ncic_min.
  Kokkos::parallel_for(
      "advance_ncm::loop25", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Loop indices
        const int i = column;
        for (k = 0; k < gr_nzt; ++k) {
          if ((Nc_in_cloud(i, k) < Ncic_min(i, k))) {
            Nc_in_cloud(i, k) = Ncic_min(i, k);
          }
        }
      });

  // Solve for < w'N_c' > at all intermediate (momentum) grid levels, using
  // a down-gradient approximation:  < w'N_c' > = - K * d< N_c >/dz.
  // A Crank-Nicholson time-stepping scheme is used for this variable.
  // This is the portion of the calculation using < N_c > from timestep t+1.
  Kokkos::parallel_for(
      "advance_ncm::loop26", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k; // Loop indices
        const int i = column;
        for (k = 0; k < gr_nzm; ++k) {
          K_Nc_nu_hm(i, k) = (K_Nc(i, k) + nu_vert_res_dep_nu_hm(i));
        }
      });

  calc_xpwp(gr.nzm, gr.nzt, ngrdcol, gr, K_Nc_nu_hm, Ncm, xpwp);

  Kokkos::parallel_for(
      "advance_ncm::array27", ((wpNcp.extent_int(0) * gr.nzm - 2 * wpNcp.extent_int(0))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0               = ((cell % wpNcp.extent_int(0)));
        const int index_1               = ((cell / wpNcp.extent_int(0)));
        wpNcp(index_0, ((index_1 + 1))) = ((-one_half) * xpwp(index_0, ((index_1 + 1))));
      });

  // A zero-flux boundary condition is used for N_c.
  Kokkos::parallel_for(
      "advance_ncm::array28", wpNcp.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0 = cell;
        wpNcp(index_0, 0) = zero;
      });
  Kokkos::parallel_for(
      "advance_ncm::array29", wpNcp.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0              = cell;
        wpNcp(index_0, ((gr_nzm - 1))) = zero;
      });

  // Statistics for all covariances involving N_c:  < w'N_c' >.

  if (stats.l_sample) {
    stats_update(std::string("wpNcp"), wpNcp, stats);
  }

  if (stats.l_sample) {
    Real2D array_actual_30 = workspace_scope.view<Real2D>("advance_ncm::array_actual_30",
                                                          Ncm.extent_int(0), Ncm.extent_int(1));
    Kokkos::parallel_for(
        "advance_ncm::array31", ((Ncm.extent_int(0) * Ncm.extent_int(1))), KOKKOS_LAMBDA(int cell) {
          const int index_0                 = ((cell % Ncm.extent_int(0)));
          const int index_1                 = ((cell / Ncm.extent_int(0)));
          array_actual_30(index_0, index_1) = (Ncm(index_0, index_1) / dt);
        });
    stats_finalize_budget(std::string("Ncm_bt"), array_actual_30, stats);
  }

  return;
}
inline void microphys_solve(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                            const std::string& solve_type, bool l_sed, Real3D lhs_ta, Real3D lhs_ma,
                            Real3D sed_turb_lhs, Real3D sed_diff_lhs, Real2D cloud_frac,
                            int tridiag_solve_method, stats_type& stats, Real3D lhs, Real2D rhs,
                            Real2D hmm, err_info_type& err_info,
                            const parameters_microphys& microphys) {
  // Description:
  // Solve the tridiagonal system for hydrometeor variable.
  // References:
  //  None
  //---------------------------------------------------------------------------
  //------------------------ Input Variables ------------------------
  //------------------------ Input/Output Variables ------------------------
  //------------------------ Local Variables ------------------------
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope          = workspace.scope();
  const auto l_in_cloud_Nc_diff = microphys.l_in_cloud_Nc_diff;
  const bool l_sample           = stats.l_sample;
  const bool l_solve_Ncm        = solve_type == "Ncm";

  std::string name_ma;
  std::string name_ta;
  std::string name_sd;
  int name_sd_length = 0;
  std::string name_ts;
  int name_ts_length = 0;
  Real2D ma_term     = workspace_scope.view<Real2D>("microphys_solve::ma_term", ngrdcol, gr.nzt);
  Real2D sd_term     = workspace_scope.view<Real2D>("microphys_solve::sd_term", ngrdcol, gr.nzt);
  Real2D ta_term     = workspace_scope.view<Real2D>("microphys_solve::ta_term", ngrdcol, gr.nzt);
  Real2D ts_term     = workspace_scope.view<Real2D>("microphys_solve::ts_term", ngrdcol, gr.nzt);
  const auto gr_nzt = gr.nzt;

  //---------------------------- Begin Code ----------------------------

  // Set up the stats variable names for the selected hydrometeor equation.
  // Initialize names to avoid compiler warnings when a branch is unused.
  name_ma        = std::string("");
  name_ta        = std::string("");
  name_sd        = std::string("");
  name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  name_ts        = std::string("");
  name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  if (solve_type == std::string("rrm")) {
    name_ma        = std::string("rrm_ma");
    name_ta        = std::string("rrm_ta");
    name_sd        = std::string("rrm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("rrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("Nrm")) {
    name_ma        = std::string("Nrm_ma");
    name_ta        = std::string("Nrm_ta");
    name_sd        = std::string("Nrm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("Nrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("rim")) {
    name_ma        = std::string("rim_ma");
    name_ta        = std::string("rim_ta");
    name_sd        = std::string("rim_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("rsm")) {
    name_ma        = std::string("rsm_ma");
    name_ta        = std::string("rsm_ta");
    name_sd        = std::string("rsm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("rgm")) {
    name_ma        = std::string("rgm_ma");
    name_ta        = std::string("rgm_ta");
    name_sd        = std::string("rgm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Ncm")) {
    name_ma = std::string("Ncm_ma");
    name_ta = std::string("Ncm_ta");
  } else if (solve_type == std::string("Nim")) {
    name_ma        = std::string("Nim_ma");
    name_ta        = std::string("Nim_ta");
    name_sd        = std::string("Nim_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Nsm")) {
    name_ma        = std::string("Nsm_ma");
    name_ta        = std::string("Nsm_ta");
    name_sd        = std::string("Nsm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Ngm")) {
    name_ma        = std::string("Ngm_ma");
    name_ta        = std::string("Ngm_ta");
    name_sd        = std::string("Ngm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else {
    name_ma        = std::string("");
    name_ta        = std::string("");
    name_sd        = std::string("");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  }

  // Solve system using a tridiag_solve.
  { // Native operators are band first; the existing solver is column first.
    Real3D lhs_solver =
        workspace_scope.view<Real3D>("microphys_solve::lhs_solver", ngrdcol, 3, gr.nzt);
    Real2D soln_solver =
        workspace_scope.view<Real2D>("microphys_solve::soln_solver", ngrdcol, gr.nzt);
    const int nzt_solver = gr.nzt;
    Kokkos::parallel_for(
        "microphys_solve::to_solver", ((ngrdcol * nzt_solver)), KOKKOS_LAMBDA(int cell) {
          const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
          for (int b = 0; b < 3; ++b)
            lhs_solver(i, b, k) = lhs(b, i, k);
        });
    tridiag_solve(workspace, solve_type, tridiag_solve_method, ngrdcol, gr.nzt, false, lhs_solver,
                  rhs, err_info, soln_solver);
    Kokkos::parallel_for(
        "microphys_solve::from_solver", ((ngrdcol * nzt_solver)), KOKKOS_LAMBDA(int cell) {
          const int i = ((cell % ngrdcol)), k = ((cell / ngrdcol));
          for (int b = 0; b < 3; ++b)
            lhs(b, i, k) = lhs_solver(i, b, k);
          hmm(i, k) = soln_solver(i, k);
        });
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << "\n";
      std::cerr << std::string("calling tridiag_solve in microphys_solve") << "\n";
      return;
    }
  }

  Kokkos::parallel_for(
      "microphys_solve::array1", ((ma_term.extent_int(0) * ma_term.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0         = ((cell % ma_term.extent_int(0)));
        const int index_1         = ((cell / ma_term.extent_int(0)));
        ma_term(index_0, index_1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "microphys_solve::array2", ((sd_term.extent_int(0) * sd_term.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0         = ((cell % sd_term.extent_int(0)));
        const int index_1         = ((cell / sd_term.extent_int(0)));
        sd_term(index_0, index_1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "microphys_solve::array3", ((ta_term.extent_int(0) * ta_term.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0         = ((cell % ta_term.extent_int(0)));
        const int index_1         = ((cell / ta_term.extent_int(0)));
        ta_term(index_0, index_1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "microphys_solve::array4", ((ts_term.extent_int(0) * ts_term.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0         = ((cell % ts_term.extent_int(0)));
        const int index_1         = ((cell / ts_term.extent_int(0)));
        ts_term(index_0, index_1) = real_t(0.0);
      });

  // Statistics: implicit contributions to hydrometeor hmm.

  if (stats.l_sample) {
    Kokkos::parallel_for(
        "microphys_solve::loop5", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k;   // Loop indices
          int km1; // Array indices
          int kp1; // Array indices
          const int i = column;
          for (k = 0; k < gr_nzt; ++k) {

            km1 = (k > 1 ? ((k - 1)) : 0);
            kp1 = (k < (gr_nzt - 2) ? ((k + 1)) : (((gr_nzt - 1))));

            if ((l_solve_Ncm && l_in_cloud_Nc_diff)) {
              ma_term(i, k) = (((((-lhs_ma(2, i, k)) * hmm(i, km1)) *
                                 ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                        : (cloud_frac_min))) -
                                ((lhs_ma(1, i, k) * hmm(i, k)) *
                                 ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                        : (cloud_frac_min)))) -
                               ((lhs_ma(0, i, k) * hmm(i, kp1)) *
                                ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                       : (cloud_frac_min))));
              ta_term(i, k) = (((((-lhs_ta(2, i, k)) * hmm(i, km1)) *
                                 ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                        : (cloud_frac_min))) -
                                ((lhs_ta(1, i, k) * hmm(i, k)) *
                                 ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                        : (cloud_frac_min)))) -
                               ((lhs_ta(0, i, k) * hmm(i, kp1)) *
                                ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                       : (cloud_frac_min))));
            } else {
              ma_term(i, k) =
                  ((((-lhs_ma(2, i, k)) * hmm(i, km1)) - (lhs_ma(1, i, k) * hmm(i, k))) -
                   (lhs_ma(0, i, k) * hmm(i, kp1)));
              ta_term(i, k) =
                  ((((-lhs_ta(2, i, k)) * hmm(i, km1)) - (lhs_ta(1, i, k) * hmm(i, k))) -
                   (lhs_ta(0, i, k) * hmm(i, kp1)));
            }

            if ((l_sed && (name_sd_length > 0))) {
              sd_term(i, k) = ((((-sed_diff_lhs(2, i, k)) * hmm(i, km1)) -
                                (sed_diff_lhs(1, i, k) * hmm(i, k))) -
                               (sed_diff_lhs(0, i, k) * hmm(i, kp1)));
            }

            if ((l_sed && (name_ts_length > 0))) {
              ts_term(i, k) = ((((-sed_turb_lhs(2, i, k)) * hmm(i, km1)) -
                                (sed_turb_lhs(1, i, k) * hmm(i, k))) -
                               (sed_turb_lhs(0, i, k) * hmm(i, kp1)));
            }
          }
        });

    if ((static_cast<int>(microphys_detail::trim(name_ma).size()) > 0)) {
      stats_update(name_ma, ma_term, stats);
    }
    if (((static_cast<int>(microphys_detail::trim(name_sd).size()) > 0) && l_sed)) {
      stats_update(name_sd, sd_term, stats);
    }
    if (((static_cast<int>(microphys_detail::trim(name_ts).size()) > 0) && l_sed)) {
      // Finalize implicit contributions.
      stats_finalize_budget(name_ts, ts_term, stats);
    }
    if ((static_cast<int>(microphys_detail::trim(name_ta).size()) > 0)) {
      // For Ncm, diffusion in this solve is formed using in-cloud values.
      // The cloud-fraction scaling is accounted for in the sampled tendency
      // so the turbulent-advection budget remains balanced.
      stats_finalize_budget(name_ta, ta_term, stats);
    }
  }

  return;
}
inline void microphys_lhs(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                          const std::string& solve_type, bool l_sed, real_t dt, Real2D K_hm,
                          Real1D nu, Real2D wm_zt, Real2D V_hm, Real2D V_hmt,
                          Real2D Vhmphmp_zt_impc, Real2D rho_ds_zm, Real2D rho_ds_zt,
                          Real2D invrs_rho_ds_zt, bool l_upwind_xm_ma, stats_type& stats,
                          Real3D lhs_ta, Real3D lhs_ma, Real3D sed_turb_lhs, Real3D sed_diff_lhs,
                          Real3D lhs, const parameters_microphys& microphys) {
  // Description:
  // Setup the matrix of implicit contributions to a term.
  // Can include the effects of sedimentation, diffusion, and advection.
  // The Morrison microphysics has an explicit sedimentation code, which is
  // handled elsewhere.
  //
  // Notes:
  // Setup for tridiagonal system and boundary conditions should be the same as
  // the original rain subroutine code.
  //-----------------------------------------------------------------------
  //------------------------------ Input Variables ------------------------------
  // Constant parameters
  // Input Variables
  // approximation rather than a centered differencing for turbulent or
  // mean advection terms. It affects rtm, thlm, sclrm, um and vm.
  //------------------------------ Output Variables ------------------------------
  //------------------------------ Local Variables ------------------------------
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope         = workspace.scope();
  const auto l_upwind_diff_sed = microphys.l_upwind_diff_sed;
  const bool l_sample          = stats.l_sample;

  constexpr int kp1_tdiag = 0;
  constexpr int k_tdiag   = 1;
  constexpr int km1_tdiag = 2;
  Real2D Vhmphmp_impc     = workspace_scope.view<Real2D>(
      "microphys_lhs::Vhmphmp_impc", ngrdcol,
      gr.nzm); // Implicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]
  Real2D Kh_zt =
      workspace_scope.view<Real2D>("microphys_lhs::Kh_zt", ngrdcol,
                                   gr.nzt); // Eddy diffusivity coefficient, thermo levels [m2/s]
  Real2D invrs_rho_ds_zt_col =
      workspace_scope.view<Real2D>("microphys_lhs::invrs_rho_ds_zt_col", ngrdcol,
                                   gr.nzt); // Inv. dry, static density @ thermo. levs.    [m^3/kg]
  Real2D Kh_zm = workspace_scope.view<Real2D>(
      "microphys_lhs::Kh_zm", ngrdcol,
      gr.nzm); // Kh Eddy diffusivity on momentum grid            [m^2/s]
  Real2D rho_ds_zm_col = workspace_scope.view<Real2D>(
      "microphys_lhs::rho_ds_zm_col", ngrdcol,
      gr.nzm); // Dry, static density on momentum levels        [kg/m^3]
  Real1D nu_col = workspace_scope.view<Real1D>("microphys_lhs::nu_col", ngrdcol);
  Real2D wm_zt_col =
      workspace_scope.view<Real2D>("microphys_lhs::wm_zt_col", ngrdcol, gr.nzt); // In
  std::string name_sd;
  int name_sd_length = 0;
  std::string name_ts;
  int name_ts_length      = 0;
  const auto gr_nzt       = gr.nzt;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  //---------------------------- Begin Code ----------------------------

  name_sd        = std::string("");
  name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  name_ts        = std::string("");
  name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  if (solve_type == std::string("rrm")) {
    name_sd        = std::string("rrm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("rrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("Nrm")) {
    name_sd        = std::string("Nrm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("Nrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("rim")) {
    name_sd        = std::string("rim_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("rsm")) {
    name_sd        = std::string("rsm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("rgm")) {
    name_sd        = std::string("rgm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Nim")) {
    name_sd        = std::string("Nim_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Nsm")) {
    name_sd        = std::string("Nsm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else if (solve_type == std::string("Ngm")) {
    name_sd        = std::string("Ngm_sd");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
  } else {
    name_sd        = std::string("");
    name_sd_length = static_cast<int>(microphys_detail::trim(name_sd).size());
    name_ts        = std::string("");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  }

  // Interpolate the implicit component of < V_hm'h_m' >, a momentum-level
  // variable that is calculated on thermodynamic levels, from thermodynamic
  // levels to momentum levels.
  zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Vhmphmp_zt_impc, Vhmphmp_impc);

  // Reset LHS Matrix for current timestep.
  Kokkos::parallel_for(
      "microphys_lhs::array1", ((lhs.extent_int(0) * lhs.extent_int(1) * lhs.extent_int(2))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0              = ((cell % lhs.extent_int(0)));
        const int index_1              = (((cell / lhs.extent_int(0)) % lhs.extent_int(1)));
        const int index_2              = ((cell / (lhs.extent_int(0) * lhs.extent_int(1))));
        lhs(index_0, index_1, index_2) = zero;
      });

  // LHS turbulent advection term.
  // - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
  // Note:  a down gradient closure approximation is made for < w'h_m' >, so
  //        the turbulent advection term is solved as an eddy-diffusion
  //        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
  // A Crank-Nicholson time-stepping scheme is used for this term.
  Kokkos::parallel_for(
      "microphys_lhs::array2", ((Kh_zm.extent_int(0) * Kh_zm.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0       = ((cell % Kh_zm.extent_int(0)));
        const int index_1       = ((cell / Kh_zm.extent_int(0)));
        Kh_zm(index_0, index_1) = K_hm(index_0, index_1);
      });
  auto interpolated_4 = zm2zt_api(gr.nzm, gr.nzt, ngrdcol, gr, K_hm);
  Kokkos::parallel_for(
      "microphys_lhs::array4", ((Kh_zt.extent_int(0) * Kh_zt.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % Kh_zt.extent_int(0)));
        const int index_1 = ((cell / Kh_zt.extent_int(0)));
        Kh_zt(index_0, index_1) =
            ((interpolated_4(index_0, index_1)) > (real_t(0.)) ? (interpolated_4(index_0, index_1))
                                                               : (real_t(0.)));
      });
  Kokkos::parallel_for(
      "microphys_lhs::array5",
      ((invrs_rho_ds_zt_col.extent_int(0) * invrs_rho_ds_zt_col.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0                     = ((cell % invrs_rho_ds_zt_col.extent_int(0)));
        const int index_1                     = ((cell / invrs_rho_ds_zt_col.extent_int(0)));
        invrs_rho_ds_zt_col(index_0, index_1) = invrs_rho_ds_zt(index_0, index_1);
      });
  Kokkos::parallel_for(
      "microphys_lhs::array6", ((rho_ds_zm_col.extent_int(0) * rho_ds_zm_col.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0               = ((cell % rho_ds_zm_col.extent_int(0)));
        const int index_1               = ((cell / rho_ds_zm_col.extent_int(0)));
        rho_ds_zm_col(index_0, index_1) = rho_ds_zm(index_0, index_1);
      });
  Kokkos::parallel_for(
      "microphys_lhs::array7", nu_col.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0 = cell;
        nu_col(index_0) = nu(index_0);
      });

  diffusion_zt_lhs(workspace, gr.nzm, gr.nzt, ngrdcol, gr, Kh_zm, Kh_zt, nu_col,
                   invrs_rho_ds_zt_col, rho_ds_zm_col, lhs_ta);

  Kokkos::parallel_for(
      "microphys_lhs::array8",
      ((lhs_ta.extent_int(0) * lhs_ta.extent_int(1) * lhs_ta.extent_int(2))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % lhs_ta.extent_int(0)));
        const int index_1 = (((cell / lhs_ta.extent_int(0)) % lhs_ta.extent_int(1)));
        const int index_2 = ((cell / (lhs_ta.extent_int(0) * lhs_ta.extent_int(1))));
        lhs_ta(index_0, index_1, index_2) = (one_half * lhs_ta(index_0, index_1, index_2));
      });

  // The lower boundary condition needs to be applied here at level 1.

  // Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
  Kokkos::parallel_for(
      "microphys_lhs::loop9", ngrdcol, KOKKOS_LAMBDA(int column) {
        const int i = column;
        lhs_ta(kp1_tdiag, i, 0) =
            (((-one_half) * invrs_rho_ds_zt(i, 0)) *
             (((gr_invrs_dzt(i, 0) * (Kh_zm(i, 1) + nu(i))) * rho_ds_zm(i, 1)) *
              gr_invrs_dzm(i, 1)));

        // Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
        lhs_ta(k_tdiag, i, 0) = (((+one_half) * invrs_rho_ds_zt(i, 0)) *
                                 (((gr_invrs_dzt(i, 0) * (Kh_zm(i, 1) + nu(i))) * rho_ds_zm(i, 1)) *
                                  gr_invrs_dzm(i, 1)));

        // Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
        lhs_ta(km1_tdiag, i, 0) = zero;
      });

  Kokkos::parallel_for(
      "microphys_lhs::array10", ((wm_zt_col.extent_int(0) * wm_zt_col.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0           = ((cell % wm_zt_col.extent_int(0)));
        const int index_1           = ((cell / wm_zt_col.extent_int(0)));
        wm_zt_col(index_0, index_1) = wm_zt(index_0, index_1);
      });

  // LHS mean advection term.
  term_ma_zt_lhs(gr.nzm, gr.nzt, ngrdcol, wm_zt_col, gr.weights_zt2zm, gr.invrs_dzt, gr.invrs_dzm,
                 l_upwind_xm_ma, gr.grid_dir, lhs_ma);

  // Setup LHS Matrix
  Kokkos::parallel_for(
      "microphys_lhs::loop22", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k;   // Loop indices
        int km1; // Array indices
        int kp1; // Array indices
        const int i = column;
        for (k = 0; k < gr_nzt; ++k) {

          km1 = (k > 1 ? ((k - 1)) : 0);
          kp1 = (k < (gr_nzt - 2) ? ((k + 1)) : (((gr_nzt - 1))));

          // LHS time tendency.
          lhs(k_tdiag, i, k) = (lhs(k_tdiag, i, k) + (one / dt));

          for (int index_0 = 0; index_0 < lhs.extent_int(0); ++index_0) {
            lhs(index_0, i, k) = (lhs(index_0, i, k) + lhs_ta(index_0, i, k));
          }

          // LHS mean advection term.
          for (int index_0 = 0; index_0 < (km1_tdiag - kp1_tdiag + 1); ++index_0) {
            lhs(((index_0 + kp1_tdiag)), i, k) =
                (lhs(((index_0 + kp1_tdiag)), i, k) + lhs_ma(((index_0 + kp1_tdiag)), i, k));
          }

          if (l_sed) {

            // LHS mean sedimentation term.
            // Note: the Morrison microphysics has its own sedimentation code,
            // which is applied through the _mc terms to each hydrometeor species.
            // Therefore, l_sed will always be false when the Morrison microphysics
            // is enabled.  -dschanen 24 Jan 2011
            if ((!l_upwind_diff_sed)) {

              // Sedimentation (both mean and turbulent) uses centered
              // differencing.  This is the default method.
              for (int index_0 = 0; index_0 < (km1_tdiag - kp1_tdiag + 1); ++index_0) {
                lhs(((index_0 + kp1_tdiag)), i, k) =
                    (lhs(((index_0 + kp1_tdiag)), i, k) +
                     sed_centered_diff_lhs(gr, V_hm(i, kp1), V_hm(i, k), rho_ds_zm(i, kp1),
                                           rho_ds_zm(i, k), invrs_rho_ds_zt(i, k),
                                           gr_invrs_dzt(i, k), k)[index_0]);
              }

            } else {

              // Sedimentation (both mean and turbulent) uses "upwind"
              // differencing.
              for (int index_0 = 0; index_0 < (km1_tdiag - kp1_tdiag + 1); ++index_0) {
                lhs(((index_0 + kp1_tdiag)), i, k) =
                    (lhs(((index_0 + kp1_tdiag)), i, k) +
                     sed_upwind_diff_lhs(gr, V_hmt(i, k), V_hmt(i, kp1), rho_ds_zt(i, k),
                                         rho_ds_zt(i, kp1), invrs_rho_ds_zt(i, k),
                                         gr_invrs_dzm(i, kp1), k)[index_0]);
              }
            }

            // LHS turbulent sedimentation term.
            for (int index_0 = 0; index_0 < (km1_tdiag - kp1_tdiag + 1); ++index_0) {
              lhs(((index_0 + kp1_tdiag)), i, k) =
                  (lhs(((index_0 + kp1_tdiag)), i, k) +
                   term_turb_sed_lhs(gr, Vhmphmp_impc(i, kp1), Vhmphmp_impc(i, k),
                                     Vhmphmp_zt_impc(i, kp1), Vhmphmp_zt_impc(i, k),
                                     rho_ds_zm(i, kp1), rho_ds_zm(i, k), rho_ds_zt(i, kp1),
                                     rho_ds_zt(i, k), gr_invrs_dzt(i, k), gr_invrs_dzm(i, kp1),
                                     invrs_rho_ds_zt(i, k), k, l_upwind_diff_sed)[index_0]);
            }

            // l_sed
          }

          if (l_sample) {

            if ((l_sed && (name_sd_length > 0))) {

              if ((!l_upwind_diff_sed)) {
                for (int index_0 = 0; index_0 < 3; ++index_0) {
                  sed_diff_lhs(index_0, i, k) = sed_centered_diff_lhs(
                      gr, V_hm(i, kp1), V_hm(i, k), rho_ds_zm(i, kp1), rho_ds_zm(i, k),
                      invrs_rho_ds_zt(i, k), gr_invrs_dzt(i, k), k)[index_0];
                }
              } else {
                for (int index_0 = 0; index_0 < 3; ++index_0) {
                  sed_diff_lhs(index_0, i, k) = sed_upwind_diff_lhs(
                      gr, V_hmt(i, k), V_hmt(i, kp1), rho_ds_zt(i, k), rho_ds_zt(i, kp1),
                      invrs_rho_ds_zt(i, k), gr_invrs_dzm(i, kp1), k)[index_0];
                }
              }
            }

            if ((l_sed && (name_ts_length > 0))) {
              for (int index_0 = 0; index_0 < 3; ++index_0) {
                sed_turb_lhs(index_0, i, k) = term_turb_sed_lhs(
                    gr, Vhmphmp_impc(i, kp1), Vhmphmp_impc(i, k), Vhmphmp_zt_impc(i, kp1),
                    Vhmphmp_zt_impc(i, k), rho_ds_zm(i, kp1), rho_ds_zm(i, k), rho_ds_zt(i, kp1),
                    rho_ds_zt(i, k), gr_invrs_dzt(i, k), gr_invrs_dzm(i, kp1),
                    invrs_rho_ds_zt(i, k), k, l_upwind_diff_sed)[index_0];
              }
            }
          }

          if ((l_sample && l_sed)) {

            // Statistics: implicit contributions to hydrometeor hmm (new stats).
            if ((name_sd_length > 0)) {

              if ((!l_upwind_diff_sed)) {
                for (int index_0 = 0; index_0 < 3; ++index_0) {
                  sed_diff_lhs(index_0, i, k) = sed_centered_diff_lhs(
                      gr, V_hm(i, kp1), V_hm(i, k), rho_ds_zm(i, kp1), rho_ds_zm(i, k),
                      invrs_rho_ds_zt(i, k), gr_invrs_dzt(i, k), k)[index_0];
                }
              } else {
                for (int index_0 = 0; index_0 < 3; ++index_0) {
                  sed_diff_lhs(index_0, i, k) = sed_upwind_diff_lhs(
                      gr, V_hmt(i, k), V_hmt(i, kp1), rho_ds_zt(i, k), rho_ds_zt(i, kp1),
                      invrs_rho_ds_zt(i, k), gr_invrs_dzm(i, kp1), k)[index_0];
                }
              }
            }

            if ((name_ts_length > 0)) {
              for (int index_0 = 0; index_0 < 3; ++index_0) {
                sed_turb_lhs(index_0, i, k) = term_turb_sed_lhs(
                    gr, Vhmphmp_impc(i, kp1), Vhmphmp_impc(i, k), Vhmphmp_zt_impc(i, kp1),
                    Vhmphmp_zt_impc(i, k), rho_ds_zm(i, kp1), rho_ds_zm(i, k), rho_ds_zt(i, kp1),
                    rho_ds_zt(i, k), gr_invrs_dzt(i, k), gr_invrs_dzm(i, kp1),
                    invrs_rho_ds_zt(i, k), k, l_upwind_diff_sed)[index_0];
              }
            }
          }

          // Native loop: i = 1, ngrdcol

          // 1..gr%nzt
        }
      });

  return;
}
inline void microphys_rhs(ScratchWorkspace& workspace, const grid& gr, int ngrdcol,
                          const std::string& solve_type, real_t dt, bool l_sed, Real2D hmm,
                          Real2D hmm_tndcy, Real2D K_hm, Real1D nu, Real2D cloud_frac,
                          Real2D Vhmphmp_zt_expc, Real2D rho_ds_zm, Real2D rho_ds_zt,
                          Real2D invrs_rho_ds_zt, stats_type& stats, Real2D rhs,
                          const parameters_microphys& microphys) {
  // Description:
  // Compute RHS vector for a given hydrometeor.
  // This subroutine computes the explicit portion of the predictive equation
  // for a given hydrometeor.
  // References:
  //-----------------------------------------------------------------------
  // Constant parameters
  //-------------------------- Input Variables --------------------------
  //-------------------------- Output Variable --------------------------
  //-------------------------- Local Variables --------------------------
  // Appended state replaces the implicit Native parameters_microphys module.
  auto workspace_scope          = workspace.scope();
  const auto l_upwind_diff_sed  = microphys.l_upwind_diff_sed;
  const auto l_in_cloud_Nc_diff = microphys.l_in_cloud_Nc_diff;
  const bool l_sample           = stats.l_sample;
  const bool l_solve_Ncm        = solve_type == "Ncm";

  constexpr int kp1_tdiag = 0;
  constexpr int k_tdiag   = 1;
  constexpr int km1_tdiag = 2;
  Real2D Vhmphmp_expc     = workspace_scope.view<Real2D>(
      "microphys_rhs::Vhmphmp_expc", ngrdcol,
      gr.nzm); // Explicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]
  Real2D Kh_zt =
      workspace_scope.view<Real2D>("microphys_rhs::Kh_zt", ngrdcol,
                                   gr.nzt); // Eddy diffusivity coefficient, thermo levels [m2/s]
  Real2D invrs_rho_ds_zt_col =
      workspace_scope.view<Real2D>("microphys_rhs::invrs_rho_ds_zt_col", ngrdcol,
                                   gr.nzt); // Inv. dry, static density @ thermo. levs.    [m^3/kg]
  Real2D Kh_zm = workspace_scope.view<Real2D>(
      "microphys_rhs::Kh_zm", ngrdcol,
      gr.nzm); // Kh Eddy diffusivity on momentum grid            [m^2/s]
  Real2D rho_ds_zm_col = workspace_scope.view<Real2D>(
      "microphys_rhs::rho_ds_zm_col", ngrdcol,
      gr.nzm); // Dry, static density on momentum levels        [kg/m^3]
  Real1D nu_col = workspace_scope.view<Real1D>("microphys_rhs::nu_col", ngrdcol);
  Real3D lhs_ta = workspace_scope.view<Real3D>(
      "microphys_rhs::lhs_ta", 3, ngrdcol,
      gr.nzt); // LHS corresponding to contribution from turbulent adv.  [1/s]
  std::string name_ta;
  std::string name_ts;
  int name_ts_length = 0;
  Real2D ta_rhs           = workspace_scope.view<Real2D>("microphys_rhs::ta_rhs", ngrdcol, gr.nzt);
  Real2D ts_rhs           = workspace_scope.view<Real2D>("microphys_rhs::ts_rhs", ngrdcol, gr.nzt);
  const auto gr_nzt       = gr.nzt;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  //-------------------------- Begin Code --------------------------

  name_ta        = std::string("");
  name_ts        = std::string("");
  name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  if (solve_type == std::string("rrm")) {
    name_ta        = std::string("rrm_ta");
    name_ts        = std::string("rrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("Nrm")) {
    name_ta        = std::string("Nrm_ta");
    name_ts        = std::string("Nrm_ts");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  } else if (solve_type == std::string("rim")) {
    name_ta = std::string("rim_ta");
  } else if (solve_type == std::string("rsm")) {
    name_ta = std::string("rsm_ta");
  } else if (solve_type == std::string("rgm")) {
    name_ta = std::string("rgm_ta");
  } else if (solve_type == std::string("Ncm")) {
    name_ta = std::string("Ncm_ta");
  } else if (solve_type == std::string("Nim")) {
    name_ta = std::string("Nim_ta");
  } else if (solve_type == std::string("Nsm")) {
    name_ta = std::string("Nsm_ta");
  } else if (solve_type == std::string("Ngm")) {
    name_ta = std::string("Ngm_ta");
  } else {
    name_ta        = std::string("");
    name_ts        = std::string("");
    name_ts_length = static_cast<int>(microphys_detail::trim(name_ts).size());
  }

  // Interpolate the explicit component of < V_hm'h_m' >, a momentum-level
  // variable that is calculated on thermodynamic levels, from thermodynamic
  // levels to momentum levels.
  zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Vhmphmp_zt_expc, Vhmphmp_expc);

  // Initialize right-hand side vector to 0.
  const int name_ta_length = static_cast<int>(microphys_detail::trim(name_ta).size());
  Kokkos::parallel_for(
      "microphys_rhs::array1", ((rhs.extent_int(0) * rhs.extent_int(1))), KOKKOS_LAMBDA(int cell) {
        const int index_0     = ((cell % rhs.extent_int(0)));
        const int index_1     = ((cell / rhs.extent_int(0)));
        rhs(index_0, index_1) = zero;
      });

  // LHS turbulent advection term.
  // - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
  // Note:  a down gradient closure approximation is made for < w'h_m' >, so
  //        the turbulent advection term is solved as an eddy-diffusion
  //        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
  // A Crank-Nicholson time-stepping scheme is used for this term.
  Kokkos::parallel_for(
      "microphys_rhs::array2", ((Kh_zm.extent_int(0) * Kh_zm.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0       = ((cell % Kh_zm.extent_int(0)));
        const int index_1       = ((cell / Kh_zm.extent_int(0)));
        Kh_zm(index_0, index_1) = K_hm(index_0, index_1);
      });
  auto interpolated_4 = zm2zt_api(gr.nzm, gr.nzt, ngrdcol, gr, K_hm);
  Kokkos::parallel_for(
      "microphys_rhs::array4", ((Kh_zt.extent_int(0) * Kh_zt.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % Kh_zt.extent_int(0)));
        const int index_1 = ((cell / Kh_zt.extent_int(0)));
        Kh_zt(index_0, index_1) =
            ((interpolated_4(index_0, index_1)) > (zero) ? (interpolated_4(index_0, index_1))
                                                         : (zero));
      });
  Kokkos::parallel_for(
      "microphys_rhs::array5",
      ((invrs_rho_ds_zt_col.extent_int(0) * invrs_rho_ds_zt_col.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0                     = ((cell % invrs_rho_ds_zt_col.extent_int(0)));
        const int index_1                     = ((cell / invrs_rho_ds_zt_col.extent_int(0)));
        invrs_rho_ds_zt_col(index_0, index_1) = invrs_rho_ds_zt(index_0, index_1);
      });
  Kokkos::parallel_for(
      "microphys_rhs::array6", ((rho_ds_zm_col.extent_int(0) * rho_ds_zm_col.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0               = ((cell % rho_ds_zm_col.extent_int(0)));
        const int index_1               = ((cell / rho_ds_zm_col.extent_int(0)));
        rho_ds_zm_col(index_0, index_1) = rho_ds_zm(index_0, index_1);
      });
  Kokkos::parallel_for(
      "microphys_rhs::array7", nu_col.extent_int(0), KOKKOS_LAMBDA(int cell) {
        const int index_0 = cell;
        nu_col(index_0) = nu(index_0);
      });

  diffusion_zt_lhs(workspace, gr.nzm, gr.nzt, ngrdcol, gr, Kh_zm, Kh_zt, nu_col,
                   invrs_rho_ds_zt_col, rho_ds_zm_col, lhs_ta);

  Kokkos::parallel_for(
      "microphys_rhs::array8",
      ((lhs_ta.extent_int(0) * lhs_ta.extent_int(1) * lhs_ta.extent_int(2))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0 = ((cell % lhs_ta.extent_int(0)));
        const int index_1 = (((cell / lhs_ta.extent_int(0)) % lhs_ta.extent_int(1)));
        const int index_2 = ((cell / (lhs_ta.extent_int(0) * lhs_ta.extent_int(1))));
        lhs_ta(index_0, index_1, index_2) = (one_half * lhs_ta(index_0, index_1, index_2));
      });

  // The lower boundary condition needs to be applied here at level 1.

  // Thermodynamic superdiagonal: [ x xm(k+1,<t+1>) ]
  Kokkos::parallel_for(
      "microphys_rhs::loop9", ngrdcol, KOKKOS_LAMBDA(int column) {
        const int i = column;
        lhs_ta(kp1_tdiag, i, 0) =
            (((-one_half) * invrs_rho_ds_zt(i, 0)) *
             (((gr_invrs_dzt(i, 0) * (Kh_zm(i, 1) + nu(i))) * rho_ds_zm(i, 1)) *
              gr_invrs_dzm(i, 1)));

        // Thermodynamic main diagonal: [ x xm(k,<t+1>) ]
        lhs_ta(k_tdiag, i, 0) = (((+one_half) * invrs_rho_ds_zt(i, 0)) *
                                 (((gr_invrs_dzt(i, 0) * (Kh_zm(i, 1) + nu(i))) * rho_ds_zm(i, 1)) *
                                  gr_invrs_dzm(i, 1)));

        // Thermodynamic subdiagonal: [ x xm(k-1,<t+1>) ]
        lhs_ta(km1_tdiag, i, 0) = zero;
      });

  // Hydrometeor right-hand side (explicit portion of the code).
  Kokkos::parallel_for(
      "microphys_rhs::array10", ((ta_rhs.extent_int(0) * ta_rhs.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0        = ((cell % ta_rhs.extent_int(0)));
        const int index_1        = ((cell / ta_rhs.extent_int(0)));
        ta_rhs(index_0, index_1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "microphys_rhs::array11", ((ts_rhs.extent_int(0) * ts_rhs.extent_int(1))),
      KOKKOS_LAMBDA(int cell) {
        const int index_0        = ((cell % ts_rhs.extent_int(0)));
        const int index_1        = ((cell / ts_rhs.extent_int(0)));
        ts_rhs(index_0, index_1) = real_t(0.0);
      });
  Kokkos::parallel_for(
      "microphys_rhs::loop12", ngrdcol, KOKKOS_LAMBDA(int column) {
        int k;   // Loop indices
        int kp1; // Array indices
        int km1; // Array indices
        const int i = column;
        for (k = 0; k < gr_nzt; ++k) {

          km1 = (k > 1 ? ((k - 1)) : 0);
          kp1 = (k < (gr_nzt - 2) ? ((k + 1)) : (((gr_nzt - 1))));

          // RHS time tendency.
          rhs(i, k) = (hmm(i, k) / dt);

          // RHS microphysics tendency term (autoconversion, accretion, evaporation,
          // etc.).
          rhs(i, k) = (rhs(i, k) + hmm_tndcy(i, k));

          rhs(i, k) = (((rhs(i, k) - (lhs_ta(km1_tdiag, i, k) * hmm(i, km1))) -
                        (lhs_ta(k_tdiag, i, k) * hmm(i, k))) -
                       (lhs_ta(kp1_tdiag, i, k) * hmm(i, kp1)));

          // RHS turbulent sedimentation term.
          if (l_sed) {
            rhs(i, k) = (rhs(i, k) + term_turb_sed_rhs(
                                         gr, Vhmphmp_expc(i, kp1), Vhmphmp_expc(i, k),
                                         Vhmphmp_zt_expc(i, kp1), Vhmphmp_zt_expc(i, k),
                                         rho_ds_zm(i, kp1), rho_ds_zm(i, k), rho_ds_zt(i, kp1),
                                         rho_ds_zt(i, k), gr_invrs_dzt(i, k), gr_invrs_dzm(i, kp1),
                                         invrs_rho_ds_zt(i, k), k, l_upwind_diff_sed));
          }

          if (l_sample) {
            if ((name_ta_length > 0)) {
              if ((l_solve_Ncm && l_in_cloud_Nc_diff)) {
                ta_rhs(i, k) = ((((lhs_ta(2, i, k) * hmm(i, km1)) *
                                  ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                         : (cloud_frac_min))) +
                                 ((lhs_ta(1, i, k) * hmm(i, k)) *
                                  ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                         : (cloud_frac_min)))) +
                                ((lhs_ta(0, i, k) * hmm(i, kp1)) *
                                 ((cloud_frac(i, k)) > (cloud_frac_min) ? (cloud_frac(i, k))
                                                                        : (cloud_frac_min))));
              } else {
                ta_rhs(i, k) = (((lhs_ta(2, i, k) * hmm(i, km1)) + (lhs_ta(1, i, k) * hmm(i, k))) +
                                (lhs_ta(0, i, k) * hmm(i, kp1)));
              }
            }
            if ((l_sed && (name_ts_length > 0))) {
              ts_rhs(i, k) = (-term_turb_sed_rhs(
                  gr, Vhmphmp_expc(i, kp1), Vhmphmp_expc(i, k), Vhmphmp_zt_expc(i, kp1),
                  Vhmphmp_zt_expc(i, k), rho_ds_zm(i, kp1), rho_ds_zm(i, k), rho_ds_zt(i, kp1),
                  rho_ds_zt(i, k), gr_invrs_dzt(i, k), gr_invrs_dzm(i, kp1), invrs_rho_ds_zt(i, k),
                  k, l_upwind_diff_sed));
              // term_turb_sed_rhs is subtracted on the model RHS, so reverse the sign
              // here when sampling the explicit turbulent-sedimentation budget term.
            }
          }
          // Native loop: i = 1, ngrdcol

          // Native loop: k = 1, gr%nzt, 1
        }
      });

  if (stats.l_sample) {
    if ((static_cast<int>(microphys_detail::trim(name_ta).size()) > 0)) {
      // Statistics: explicit contributions for the hydrometeor.
      // hmm term ta has both implicit and explicit components; call
      // stats_begin_budget for the explicit turbulent advection component.
      stats_begin_budget(name_ta, ta_rhs, stats);
    }
    if (((static_cast<int>(microphys_detail::trim(name_ts).size()) > 0) && l_sed)) {
      // hmm term ts has both implicit and explicit components; call
      // stats_begin_budget for the explicit turbulent sedimentation component.
      stats_begin_budget(name_ts, ts_rhs, stats);
    }
  }

  return;
}
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
sed_centered_diff_lhs(const grid& gr, real_t V_hmp1, real_t V_hm, real_t rho_ds_zmp1,
                      real_t rho_ds_zm, real_t invrs_rho_ds_zt, real_t invrs_dzt, int level) {
  // Description:
  // Mean sedimentation of a hydrometeor:  implicit portion of the code, using
  // the centered difference approximation to the vertical derivative.
  //
  // The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
  // stands for the sedimentation velocity of the aforementioned hydrometeor.
  //
  // The d(hm)/dt equation contains a sedimentation term:
  //
  // - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
  //
  // The variables hm and V_hm in the sedimentation term are divided into mean
  // and turbulent components, and the term is averaged, resulting in:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // The mean sedimentation term in the d<hm>/dt equation is:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of the term is changed to
  //        a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d<hm>/dt equation.
  //
  // This term is discretized as follows when using the centered-difference
  // approximation:
  //
  // The values of <hm> are found on the thermodynamic levels, while the values
  // of <V_hm> are found on the momentum levels.  Additionally, the values of
  // rho_ds_zm are found on the momentum levels, and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  The variable <hm>
  // is interpolated to the intermediate momentum levels.  At the intermediate
  // momentum levels, the interpolated values of <hm> are multiplied by the
  // values of <V_hm> and the values of rho_ds_zm.  Then, the derivative of
  // (rho_ds*<V_hm>*<hm>) is taken over the central thermodynamic level, where
  // it is multiplied by invrs_rho_ds_zt.
  //
  // -----hmp1------------------------------------------------ t(k+1)
  //
  // =============hm(interp)=====V_hm=====rho_ds_zm=========== m(k+1)
  //
  // -----hm--------invrs_rho_ds_zt----d(rho_ds*V_hm*hm)/dz--- t(k)
  //
  // =============hm(interp)=====V_hmm1===rho_ds_zmm1========= m(k)
  //
  // -----hmm1------------------------------------------------ t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1),
  // respectively.  The letter "t" is used for thermodynamic levels and the
  // letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  //
  // Conservation Properties:
  //
  // When a hydrometeor is sedimented to the ground (or out the lower boundary
  // of the model), it is removed from the atmosphere (or from the model
  // domain).  Thus, the quantity of the hydrometeor over the entire vertical
  // domain should not be conserved due to the process of sedimentation.  Thus,
  // not all of the column totals in the left-hand side matrix should be equal
  // to 0. Instead, the sum of all the column totals should equal the flux of
  // <hm> out the bottom (zm(1) level) of the domain,
  // -rho_ds_zm(1) * V_hm(1) * hm(1), where the value of the hydrometeor at
  // the surface (zm level 1) is set equal to the value of the hydrometeor at
  // zt level 1, which is hm(1). Furthermore, most of the individual column
  // totals should sum to 0, but the 1st and 2nd (from the left) columns should
  // combine to sum to the flux out the bottom of the domain.
  //
  // To see that this modified conservation law is satisfied, compute the
  // sedimentation of hm and integrate vertically.  In discretized matrix
  // notation (where "i" stands for the matrix column and "j" stands for the
  // matrix row):
  //
  // - rho_ds_zm(1) * V_hm(1) * hm(1)
  // = Sum_j Sum_i
  //   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzt )_i
  //   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
  //
  // The left-hand side matrix,
  // ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
  // partially written below.  The sum over i in the above equation removes
  // invrs_rho_ds_zt and invrs_dzt everywhere from the matrix below.  The sum
  // over j leaves the column totals and the flux at zm(1) that are desired.
  //
  // Left-hand side matrix contributions from the sedimentation term (only);
  // first four vertical levels:
  //
  //     -------------------------------------------------------------------->
  //k=1 |   +invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)            0
  //    |    *invrs_dzt(k)        *invrs_dzt(k)
  //    |    *[ rho_ds_zm(k+1)    *rho_ds_zm(k+1)
  //    |       *V_hm(k+1)*B(k)   *V_hm(k+1)*A(k)
  //    |      -rho_ds_zm(k)
  //    |       *V_hm(k) ]
  //    |
  //k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
  //    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
  //    |    *rho_ds_zm(k)        *[ rho_ds_zm(k+1)      *rho_ds_zm(k+1)
  //    |    *V_hm(k)*D(k)           *V_hm(k+1)*B(k)     *V_hm(k+1)*A(k)
  //    |                           -rho_ds_zm(k)
  //    |                            *V_hm(k)*C(k) ]
  //    |
  //k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
  //    |                         *invrs_dzt(k)          *invrs_dzt(k)
  //    |                         *rho_ds_zm(k)          *[ rho_ds_zm(k+1)
  //    |                         *V_hm(k)*D(k)             *V_hm(k+1)*B(k)
  //    |                                                  -rho_ds_zm(k)
  //    |                                                   *V_hm(k)*C(k) ]
  //    |
  //k=4 |           0                     0             -invrs_rho_ds_zt(k)
  //    |                                                *invrs_dzt(k)
  //    |                                                *rho_ds_zm(k)
  //    |                                                *V_hm(k)*D(k)
  //    |
  //   \ /
  //
  // The variables A(k), B(k), C(k), and D(k) are weights of interpolation
  // around the central thermodynamic level (k), such that:
  //
  // A(k) = ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ),
  // B(k) = 1 - [ ( zm(k+1) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
  //      = 1 - A(k);
  // C(k) = ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
  // D(k) = 1 - [ ( zm(k) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
  //      = 1 - C(k).
  //
  // Furthermore, for all intermediate thermodynamic grid levels (as long as
  // k /= gr%nz and k /= 1), the four weighting factors have the following
  // relationships:  A(k) = C(k+1) and B(k) = D(k+1).
  //
  // Note:  The superdiagonal term from level 3 and both the main diagonal
  //        and superdiagonal terms from level 4 are not shown on this
  //        diagram.
  // References:
  // None
  // Notes:
  //   Both COAMPS Microphysics and Brian Griffin's implementation use
  //   Khairoutdinov and Kogan (2000) for the calculation of rain
  //   mixing ratio and rain droplet number concentration sedimentation
  //   velocities, but COAMPS has only the local parameterization.
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable
  // Local Variables

  constexpr int kp1_tdiag = 0;
  constexpr int k_tdiag   = 1;
  constexpr int km1_tdiag = 2;
  constexpr int t_above   = 0;
  constexpr int t_below   = 1;
  Kokkos::Array<real_t, 3> lhs;
  int mkp1; // Momentum level directly above central thermodynamic level.
  int mk;   // Momentum level directly below central thermodynamic level.
  int i;    // Loop indices
  const auto gr_nzt           = gr.nzt;
  const auto gr_weights_zt2zm = gr.weights_zt2zm;

  // ---- Begin Code ----

  // Some procedures expect arrays that are dimension(ngrdcol,nzt), so the aruments
  // to be passed to those procedures have a dummy dimension hardcoded to 1.
  // For visual clarity and consistency with other routines, we use i to index
  // those dummy dimensions, but we always want i=1.
  i = 0;

  // Momentum level (k+1) is between thermodynamic level (k+1)
  // and thermodynamic level (k).
  mkp1 = ((level + 1));
  // Momentum level (k) is between thermodynamic level (k)
  // and thermodynamic level (k-1).
  mk = level;

  if ((level == 0)) {

    // Native loop: k = 1 (bottom level); lower boundary level.
    // Special discretization where the value of the hydrometeor at the lower
    // boundary or surface (momentum level 1) is implicitly set equal to the
    // value of the hydrometeor at thermodynamic level 1.

    // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
    lhs[kp1_tdiag] = (((((+invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zmp1) * V_hmp1) *
                      gr_weights_zt2zm(i, mkp1, t_above));

    // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
    lhs[k_tdiag] =
        (((+invrs_rho_ds_zt) * invrs_dzt) *
         (((rho_ds_zmp1 * V_hmp1) * gr_weights_zt2zm(i, mkp1, t_below)) - (rho_ds_zm * V_hm)));

    // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
    lhs[km1_tdiag] = zero;

  } else {
    if (((level > 0) && (level < gr_nzt - 1))) {

      // Most of the interior model; normal conditions.

      // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs[kp1_tdiag] = (((((+invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zmp1) * V_hmp1) *
                        gr_weights_zt2zm(i, mkp1, t_above));

      // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs[k_tdiag] = (((+invrs_rho_ds_zt) * invrs_dzt) *
                      (((rho_ds_zmp1 * V_hmp1) * gr_weights_zt2zm(i, mkp1, t_below)) -
                       ((rho_ds_zm * V_hm) * gr_weights_zt2zm(i, mk, t_above))));

      // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs[km1_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * V_hm) *
                        gr_weights_zt2zm(i, mk, t_below));

    } else {
      if ((level == gr_nzt - 1)) {

        // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
        // Special discretization where the value of the hydrometeor at the upper
        // boundary (momentum level gr%nzm) is implicitly set equal to 0. The
        // process of sedimentation can only decrease the value of the hydrometeor
        // at level gr%nzt.

        // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
        lhs[kp1_tdiag] = zero;

        // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
        lhs[k_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * V_hm) *
                        gr_weights_zt2zm(i, mk, t_above));

        // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
        lhs[km1_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * V_hm) *
                          gr_weights_zt2zm(i, mk, t_below));
      }
    }
  }

  return lhs;
}
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
sed_upwind_diff_lhs(const grid& gr, real_t V_hmt, real_t V_hmtp1, real_t rho_ds_zt,
                    real_t rho_ds_ztp1, real_t invrs_rho_ds_zt, real_t invrs_dzmp1, int level) {
  // Description:
  // Mean sedimentation of a hydrometeor:  implicit portion of the code, using
  // the "upwind" difference approximation to the vertical derivative.
  //
  // The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
  // stands for the sedimentation velocity of the aforementioned hydrometeor.
  //
  // The d(hm)/dt equation contains a sedimentation term:
  //
  // - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
  //
  // The variables hm and V_hm in the sedimentation term are divided into mean
  // and turbulent components, and the term is averaged, resulting in:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // The mean sedimentation term in the d<hm>/dt equation is:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
  //
  // This term is solved for completely implicitly, such that:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of the term is changed to
  //        a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d<hm>/dt equation.
  //
  // This term is discretized as follows when using the upwind-difference
  // approximation:
  //
  // The values of <hm> and the values of V_hmt are found on the thermodynamic
  // levels.  Additionally, the values of rho_ds_zt and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  At the
  // thermodynamic levels, the values of <hm> are multiplied by the values of
  // V_hmt and the values of rho_ds_zt.  Then, the derivative of
  // (rho_ds*<V_hm>*<hm>) is taken between the thermodynamic level above the
  // central thermodynamic level and the central thermodynamic level.  The
  // derivative is multiplied by invrs_rho_ds_zt.
  //
  // --hmp1--V_hmtp1--rho_ds_ztp1--------------------------------------- t(k+1)
  //
  // =================================================================== m(k+1)
  //
  // --hm----V_hmt----rho_ds_zt--invrs_rho_ds_zt--d(rho_ds*V_hm*hm)/dz-- t(k)
  //
  // The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
  // zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
  //
  //
  // Conservation Properties:
  //
  // When a hydrometeor is sedimented to the ground (or out the lower boundary
  // of the model), it is removed from the atmosphere (or from the model
  // domain).  Thus, the quantity of the hydrometeor over the entire vertical
  // domain should not be conserved due to the process of sedimentation.  Thus,
  // not all of the column totals in the left-hand side matrix should be equal
  // to 0. Instead, the sum of all the column totals should equal the flux of
  // <hm> out the bottom (zm(1) level, for which the value of the hydrometeor
  // is set equal to the value of the hydrometeor at the zt(1) level) of the
  // domain, -rho_ds_zt(1) * V_hmt(1) * hm(1).  Furthermore, most of the
  // individual column totals should sum to 0, but the 2nd (from the left)
  // column should be equal to the flux out the bottom of the domain.
  //
  // To see that this modified conservation law is satisfied, compute the
  // sedimentation of hm and integrate vertically.  In discretized matrix
  // notation (where "i" stands for the matrix column and "j" stands for the
  // matrix row):
  //
  // - rho_ds_zt(1) * V_hmt(1) * hm(1)
  // = Sum_j Sum_i
  //   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzm )_i
  //   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
  //
  // The left-hand side matrix,
  // ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
  // partially written below.  The sum over i in the above equation removes
  // invrs_rho_ds_zt and invrs_dzm everywhere from the matrix below.  The sum
  // over j leaves the column totals and the flux at zt(1) that are desired.
  //
  // Left-hand side matrix contributions from the sedimentation term (only);
  // first three vertical levels:
  //
  //     -------------------------------------------------------------------->
  //k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
  //    |  *invrs_dzm(k+1)        *invrs_dzm(k+1)
  //    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
  //    |  *V_hmt(k)              *V_hmt(k+1)
  //    |
  //k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
  //    |                         *invrs_dzm(k+1)        *invrs_dzm(k+1)
  //    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
  //    |                         *V_hmt(k)              *V_hmt(k+1)
  //    |
  //k=3 |           0                     0             -invrs_rho_ds_zt(k)
  //    |                                                *invrs_dzm(k+1)
  //    |                                                *rho_ds_zt(k)
  //    |                                                *V_hmt(k)
  //   \ /
  //
  // Note:  The superdiagonal term from level 3 is not shown on this diagram.
  // References:
  // None
  // Notes:
  // Both COAMPS Microphysics and Brian Griffin's implementation use
  // Khairoutdinov and Kogan (2000) for the calculation of rain
  // mixing ratio and rain droplet number concentration sedimentation
  // velocities, but COAMPS has only the local parameterization.
  //
  // Please note that "upwind" sedimentation is only 1st-order accurate and
  // highly diffusive.
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable

  constexpr int kp1_tdiag = 0;
  constexpr int k_tdiag   = 1;
  constexpr int km1_tdiag = 2;
  Kokkos::Array<real_t, 3> lhs;
  const auto gr_nzt = gr.nzt;

  // ---- Begin Code ----

  // Sedimention is always a downward process, so we omit the upward case
  // (i.e. the V_hmt variable will always be negative).
  if ((level == gr_nzt - 1)) {

    // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
    // Special discretization where the value of the hydrometeor above the
    // upper boundary is implicitly set equal to 0. The process of
    // sedimentation can only decrease the value of the hydrometeor
    // at level gr%nzt.

    // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
    lhs[kp1_tdiag] = zero;

    // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
    lhs[k_tdiag] = ((((-invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_zt) * V_hmt);

    // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
    lhs[km1_tdiag] = zero;

  } else {

    // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
    lhs[kp1_tdiag] = ((((+invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_ztp1) * V_hmtp1);

    // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
    lhs[k_tdiag] = ((((-invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_zt) * V_hmt);

    // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
    lhs[km1_tdiag] = zero;
  }

  return lhs;
}
KOKKOS_INLINE_FUNCTION Kokkos::Array<real_t, 3>
term_turb_sed_lhs(const grid& gr, real_t Vhmphmp_impcp1, real_t Vhmphmp_impc,
                  real_t Vhmphmp_zt_impcp1, real_t Vhmphmp_zt_impc, real_t rho_ds_zmp1,
                  real_t rho_ds_zm, real_t rho_ds_ztp1, real_t rho_ds_zt, real_t invrs_dzt,
                  real_t invrs_dzmp1, real_t invrs_rho_ds_zt, int level, bool l_upwind_diff_sed) {
  // Description:
  // Turbulent sedimentation of a hydrometeor:  implicit portion of the code.
  //
  // The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
  // stands for the sedimentation velocity of the aforementioned hydrometeor.
  //
  // The d(hm)/dt equation contains a sedimentation term:
  //
  // - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
  //
  // The variables hm and V_hm in the sedimentation term are divided into mean
  // and turbulent components, and the term is averaged, resulting in:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // The turbulent sedimentation term in the d<hm>/dt equation is:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
  // on < hm > in the manner:
  //
  // < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
  //
  // This term can also be solved for completely explicitly (it's original
  // form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
  // < V_hm'hm' >.  The equation becomes:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
  //
  // where the timestep index (t+1) means that the value of < hm > being used
  // is from the next timestep, which is being advanced to in solving the
  // d<hm>/dt equation.  Implicit and explicit portions of this term are
  // produced.  The implicit portion of this term is:
  //
  // - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < hm >(t+1) ) / dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of the d[ ] / dz term is
  //        changed to a "+".
  //
  // This term can be discretized using the centered-difference approximation
  // (which is preferred), or else using the "upwind"-difference approximation.
  //
  // The implicit portion of this term is discretized as follows when using
  // the centered-difference approximation:
  //
  // The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
  // thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
  // thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
  // the momentum levels, and the values of invrs_rho_ds_zt are found on the
  // thermodynamic levels.  The variables < hm > and Vhmphmp_zt_impc are both
  // interpolated to the intermediate momentum levels.  At the momentum levels,
  // the values of interpolated < hm > and interpolated Vhmphmp_zt_impc are
  // multiplied together, and their products are multiplied by the values of
  // rho_ds_zm.  The mathematical expression F is the product of these three
  // variables at momentum levels.  Then, the derivative dF/dz is taken over
  // the central thermodynamic level, where it is multiplied by
  // invrs_rho_ds_zt.  In this function, the value of F is as follows:
  //
  // F = rho_ds_zm * Vhmphmp_impc(interp) * hmm(interp).
  //
  //
  // ----hmmp1--------Vhmphmp_zt_impcp1--------------------------------- t(k+1)
  //
  // =====hmm(interp)=====Vhmphmp_impcp1(interp)=====rho_ds_zmp1======== m(k+1)
  //
  // ----hmm----------Vhmphmp_zt_impc-----invrs_rho_ds_zt-----dF/dz----- t(k)
  //
  // =====hmm(interp)=====Vhmphmp_impc(interp)=======rho_ds_zm========== m(k)
  //
  // ----hmmm1--------Vhmphmp_zt_impcm1--------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
  //
  // The implicit portion of this term is discretized as follows when using
  // the upwind-difference approximation:
  //
  // The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
  // thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
  // thermodynamic levels.  Additionally, the values of rho_ds_zt and the
  // values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
  // thermodynamic levels, the values of < hm > and Vhmphmp_zt_impc are
  // multiplied together, and their products are multiplied by the values of
  // rho_ds_zt.  The mathematical expression F is the product of these three
  // variables at thermodynamic levels.  Then, the derivative dF/dz is taken
  // between the thermodynamic level above the central thermodynamic level and
  // the central thermodynamic level.  The derivative is multiplied
  // by invrs_rho_ds_zt.  In this function, the value of F is as follows:
  //
  // F = rho_ds_zt * Vhmphmp_zt_impc * hmm.
  //
  //
  // --hmmp1---Vhmphmp_zt_impcp1---rho_ds_ztp1-------------------------- t(k+1)
  //
  // =================================================================== m(k+1)
  //
  // --hmm-----Vhmphmp_zt_impc-----rho_ds_zt---invrs_rho_ds_zt---dF/dz-- t(k)
  //
  // The vertical indices t(k+1), m(k+1), and t(k) correspond with altitudes
  // zt(k+1), zm(k+1), and zt(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).
  // References:
  //  None
  //
  // Notes:
  // Please note that "upwind" sedimentation is only 1st-order accurate and
  // highly diffusive.
  //-----------------------------------------------------------------------
  // Constant parameters
  // Input Variables
  // Return Variable
  // Local Variables

  constexpr int kp1_tdiag = 0;
  constexpr int k_tdiag   = 1;
  constexpr int km1_tdiag = 2;
  constexpr int t_above   = 0;
  constexpr int t_below   = 1;
  Kokkos::Array<real_t, 3> lhs;
  int mkp1; // Momentum level directly above central thermodynamic level.
  int mk;   // Momentum level directly below central thermodynamic level.
  int i;    // Loop indices
  const auto gr_nzt           = gr.nzt;
  const auto gr_weights_zt2zm = gr.weights_zt2zm;

  //---------------------------- Begin Code ----------------------------

  // Some procedures expect arrays that are dimension(ngrdcol,nz), so the aruments
  // to be passed to those procedures have a dummy dimension hardcoded to 1.
  // For visual clarity and consistency with other routines, we use i to index
  // those dummy dimensions, but we always want i=1.
  i = 0;

  // Momentum level (k+1) is between thermodynamic level (k+1)
  // and thermodynamic level (k).
  mkp1 = ((level + 1));
  // Momentum level (k) is between thermodynamic level (k)
  // and thermodynamic level (k-1).
  mk = level;

  // LHS (implicit component) of turbulent sedimentation term,
  // - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
  // = - (1/rho_ds)
  //     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
  // Implicit component:
  // - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < h_m > ) / dz.
  if ((!l_upwind_diff_sed)) {

    // Sedimentation (both mean and turbulent) uses centered differencing.

    if ((level == 0)) {

      // Native loop: k = 1 (bottom level); lower boundary level.
      // Special discretization where the value of the hydrometeor at the
      // lower boundary or surface (momentum level 1) is implicitly set equal
      // to the value of the hydrometeor at thermodynamic level 1.

      // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs[kp1_tdiag] = (((((+invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zmp1) * Vhmphmp_impcp1) *
                        gr_weights_zt2zm(i, mkp1, t_above));

      // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs[k_tdiag] = (((+invrs_rho_ds_zt) * invrs_dzt) *
                      (((rho_ds_zmp1 * Vhmphmp_impcp1) * gr_weights_zt2zm(i, mkp1, t_below)) -
                       (rho_ds_zm * Vhmphmp_impc)));

      // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs[km1_tdiag] = zero;

    } else {
      if (((level > 0) && (level < gr_nzt - 1))) {

        // Most of the interior model; normal conditions.

        // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
        lhs[kp1_tdiag] = (((((+invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zmp1) * Vhmphmp_impcp1) *
                          gr_weights_zt2zm(i, mkp1, t_above));

        // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
        lhs[k_tdiag] = (((+invrs_rho_ds_zt) * invrs_dzt) *
                        (((rho_ds_zmp1 * Vhmphmp_impcp1) * gr_weights_zt2zm(i, mkp1, t_below)) -
                         ((rho_ds_zm * Vhmphmp_impc) * gr_weights_zt2zm(i, mk, t_above))));

        // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
        lhs[km1_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * Vhmphmp_impc) *
                          gr_weights_zt2zm(i, mk, t_below));

      } else {
        if ((level == gr_nzt - 1)) {

          // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
          // Special discretization where the values of both the hydrometeor and
          // V_hm'hm' at the upper boundary (momentum level gr%nzm) are set equal
          // to 0.

          // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs[kp1_tdiag] = zero;

          // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs[k_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * Vhmphmp_impc) *
                          gr_weights_zt2zm(i, mk, t_above));

          // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs[km1_tdiag] = (((((-invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * Vhmphmp_impc) *
                            gr_weights_zt2zm(i, mk, t_below));

          // level
        }
      }
    }

  } else {
    // l_upwind_diff_sed

    // Sedimentation (both mean and turbulent) uses "upwind" differencing.
    if ((level < gr_nzt - 1)) {

      // Most of the interior model; normal conditions.

      // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs[kp1_tdiag] = ((((+invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_ztp1) * Vhmphmp_zt_impcp1);

      // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs[k_tdiag] = ((((-invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_zt) * Vhmphmp_zt_impc);

      // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs[km1_tdiag] = zero;

    } else {
      if ((level == gr_nzt - 1)) {

        // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
        // Special discretization where the values of both the hydrometeor and
        // V_hm'hm' at the upper boundary (momentum level gr%nzm) are set equal
        // to 0.

        // Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
        lhs[kp1_tdiag] = zero;

        // Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
        lhs[k_tdiag] = ((((-invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_zt) * Vhmphmp_zt_impc);

        // Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
        lhs[km1_tdiag] = zero;

        // level
      }
    }
  }

  return lhs;
}
KOKKOS_INLINE_FUNCTION real_t term_turb_sed_rhs(const grid& gr, real_t Vhmphmp_expcp1,
                                                real_t Vhmphmp_expc, real_t Vhmphmp_zt_expcp1,
                                                real_t Vhmphmp_zt_expc, real_t rho_ds_zmp1,
                                                real_t rho_ds_zm, real_t rho_ds_ztp1,
                                                real_t rho_ds_zt, real_t invrs_dzt,
                                                real_t invrs_dzmp1, real_t invrs_rho_ds_zt,
                                                int level, bool l_upwind_diff_sed) {
  // Description:
  // Turbulent sedimentation of a hydrometeor:  explicit portion of the code.
  //
  // The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
  // stands for the sedimentation velocity of the aforementioned hydrometeor.
  //
  // The d(hm)/dt equation contains a sedimentation term:
  //
  // - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
  //
  // The variables hm and V_hm in the sedimentation term are divided into mean
  // and turbulent components, and the term is averaged, resulting in:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // The turbulent sedimentation term in the d<hm>/dt equation is:
  //
  // - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
  //
  // This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
  // on < hm > in the manner:
  //
  // < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
  //
  // This term can also be solved for completely explicitly (it's original
  // form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
  // < V_hm'hm' >.  The equation becomes:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
  //
  // where the timestep index (t+1) means that the value of < hm > being used
  // is from the next timestep, which is being advanced to in solving the
  // d<hm>/dt equation.  Implicit and explicit portions of this term are
  // produced.  The explicit portion of this term is:
  //
  // - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
  //
  // This term can be discretized using the centered-difference approximation
  // (which is preferred), or else using the "upwind"-difference approximation.
  //
  // The explicit portion of this term is discretized as follows when using
  // the centered-difference approximation:
  //
  // The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
  // thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
  // thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
  // the momentum levels, and the values of invrs_rho_ds_zt are found on the
  // thermodynamic levels.  The variable Vhmphmp_zt_expc is interpolated to the
  // intermediate momentum levels.  At the momentum levels, the values of
  // interpolated Vhmphmp_zt_expc are multiplied by the values of rho_ds_zm.
  // Then, the derivative d(rho_ds*Vhmphmp_zt_expc)/dz is taken over the
  // central thermodynamic level, where it is multiplied by invrs_rho_ds_zt.
  //
  // ---Vhmphmp_zt_expcp1----------------------------------------------- t(k+1)
  //
  // ======Vhmphmp_expcp1(interp)=======rho_ds_zmp1===================== m(k+1)
  //
  // ---Vhmphmp_zt_expc--invrs_rho_ds_zt--d(rho_ds*Vhmphmp_zt_expc)/dz-- t(k)
  //
  // ======Vhmphmp_expc(interp)=========rho_ds_zm======================= m(k)
  //
  // ---Vhmphmp_zt_expcm1----------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).
  //
  // The explicit portion of this term is discretized as follows when using
  // the upwind-difference approximation:
  //
  // The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
  // thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
  // thermodynamic levels.  Additionally, the values of rho_ds_zt and the
  // values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
  // thermodynamic levels, the values of Vhmphmp_zt_expc are multiplied by the
  // values of rho_ds_zt.  The mathematical expression F is the product of
  // these variables at thermodynamic levels.  Then, the derivative dF/dz is
  // taken between the thermodynamic level above the central thermodynamic
  // level and the central thermodynamic level.  The derivative is multiplied
  // by invrs_rho_ds_zt.  In this function, the value of F is as follows:
  //
  // F = rho_ds_zt * Vhmphmp_zt_expc.
  //
  //
  // -----Vhmphmp_zt_expcp1---rho_ds_ztp1------------------------------- t(k+1)
  //
  // =================================================================== m(k+1)
  //
  // -----Vhmphmp_zt_expc-----rho_ds_zt----invrs_rho_ds_zt----dF/dz----- t(k)
  //
  // The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
  // zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) ).
  // References:
  //  None
  //
  // Notes:
  // Please note that "upwind" sedimentation is only 1st-order accurate and
  // highly diffusive.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  real_t rhs; // Right hand side vector
  const auto gr_nzt = gr.nzt;

  //---------------------------- Begin Code ----------------------------

  // Initialize Right-hand side term.
  rhs = zero;

  // RHS (explicit component) of turbulent sedimentation term,
  // - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
  // = - (1/rho_ds)
  //     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
  // Explicit component:  - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
  if ((!l_upwind_diff_sed)) {

    // Sedimentation (both mean and turbulent) uses centered differencing.
    if ((level < gr_nzt - 1)) {

      // Most of the interior model; normal conditions.
      rhs = (((-invrs_rho_ds_zt) * invrs_dzt) *
             ((rho_ds_zmp1 * Vhmphmp_expcp1) - (rho_ds_zm * Vhmphmp_expc)));

    } else {
      if ((level == gr_nzt - 1)) {

        // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
        // Special discretization where the value of V_hm'hm' at the upper
        // boundary (momentum level gr%nzm) is set equal to 0.
        // This results in an equation that is the same as setting
        // Vhmphmp_expcp1 to 0 as found in the equation above for most of the
        // interior model.
        rhs = ((((+invrs_rho_ds_zt) * invrs_dzt) * rho_ds_zm) * Vhmphmp_expc);
      }
    }

  } else {
    // l_upwind_diff_sed

    // Sedimentation (both mean and turbulent) uses "upwind" differencing.
    if ((level < gr_nzt - 1)) {

      // Most of the interior model; normal conditions.
      rhs = (((-invrs_rho_ds_zt) * invrs_dzmp1) *
             ((rho_ds_ztp1 * Vhmphmp_zt_expcp1) - (rho_ds_zt * Vhmphmp_zt_expc)));

    } else {
      if ((level == gr_nzt - 1)) {

        // Native loop: k = gr%nzt (top level); upper boundary level; no flux.
        // Special discretization where the value of V_hm'hm' at the upper
        // boundary (momentum level gr%nzm) is set equal to 0.
        // This results in an equation that is the same as setting
        // Vhmphmp_zt_expcp1 to 0 as found in the equation above for most of the
        // interior model.
        rhs = ((((+invrs_rho_ds_zt) * invrs_dzmp1) * rho_ds_zt) * Vhmphmp_zt_expc);
      }
    }
  }

  return rhs;
}
inline Real3D calculate_K_hm(const grid& gr, int ngrdcol, Real2D wp2, Real2D Kh_zm, Real2D Skw_zm,
                             Real2D Lscale, int hydromet_dim, Real1D hydromet_tol, Real3D hydromet,
                             Real3D hydrometp2, Real2D clubb_params,
                             bool l_use_non_local_diff_fac) {
  // Description:
  // The predictive equation for a hydrometeor, hm, contains a turbulent
  // advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'hm'> )/dz.
  //
  // The value of <w'hm'> can be calculated in many ways.  For simplicity, a
  // down-gradient approximation is used here, where:
  //
  // <w'hm'> = -K_hm * d<hm>/dz.
  //
  // The coefficient of diffusion, K_hm, is variable and depends on multiple
  // factors.  It optionally includes a non-local turbulent advection factor.
  // References:
  // CLUBB ticket 651 and CLUBB ticket 739.
  //-----------------------------------------------------------------------
  // Input Variables
  // eddy-diffusivity applied to hydrometeors.
  // Return Variable
  // Local Variables

  Real3D K_hm("calculate_K_hm::K_hm", ngrdcol, gr.nzm,
              hydromet_dim); // hm eddy diffusivity on momentum grid        [m^2/s]
  Real3D hydromet_zm("calculate_K_hm::hydromet_zm", ngrdcol, gr.nzm,
                     hydromet_dim); // Mean of hydrometeor interp. to m-levs. [units vary]
  Real2D Lscale_zm("calculate_K_hm::Lscale_zm", ngrdcol,
                   gr.nzm); // Length scale interpolated to momentum levels
  const auto gr_nzm       = gr.nzm;
  const auto gr_nzt       = gr.nzt;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  //---------------------------- Begin Code ----------------------------

  zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Lscale, Lscale_zm);
  // Each interpolation is a host API launch; species retain the Native order.
  for (int h = 0; h < hydromet_dim; ++h) {
    zt2zm_api(gr.nzm, gr.nzt, ngrdcol, gr, Kokkos::subview(hydromet, Kokkos::ALL, Kokkos::ALL, h),
              Kokkos::subview(hydromet_zm, Kokkos::ALL, Kokkos::ALL, h));
  }

  // Loop over all hydrometeors.
  Kokkos::parallel_for(
      "calculate_k_hm::loop2", ngrdcol, KOKKOS_LAMBDA(int column) {
        real_t K_gamma; // Non-local factor of diffusion (t. adv.) for hydrometeors [m^2/s]
        int k;          // Loop indices
        int km1;        // Array indices
        int h;          // Loop indices
        const int i = column;
        for (h = 0; h < hydromet_dim; ++h) {

          // Loop over all vertical levels for each hydrometeor.
          for (k = 1; k < gr_nzm - 1; ++k) {

            km1 = (k > 1 ? ((k - 1)) : 0);

            K_hm(i, k, h) = (((clubb_params(i, ic_K_hm) * Kh_zm(i, k)) *
                              (Kokkos::sqrt(hydrometp2(i, k, h)) /
                               ((hydromet_zm(i, k, h)) > (hydromet_tol(h)) ? (hydromet_zm(i, k, h))
                                                                           : (hydromet_tol(h))))) *
                             (one + Kokkos::abs(Skw_zm(i, k))));

            if (l_use_non_local_diff_fac) {
              K_gamma =
                  (one -
                   (clubb_params(i, ic_K_hmb) *
                    ((((Lscale_zm(i, k)) > (real_t(0.0)) ? (Lscale_zm(i, k)) : (real_t(0.0))) /
                      ((hydromet_zm(i, k, h)) > (hydromet_tol(h)) ? (hydromet_zm(i, k, h))
                                                                  : (hydromet_tol(h)))) *
                     (gr_invrs_dzm(i, k) * (hydromet(i, k, h) - hydromet(i, km1, h))))));

              K_hm(i, k, h) = (K_hm(i, k, h) * ((K_gamma) > (clubb_params(i, iK_hm_min_coef))
                                                    ? (K_gamma)
                                                    : (clubb_params(i, iK_hm_min_coef))));
            }

            if ((Kokkos::abs((gr_invrs_dzm(i, k) * (hydromet(i, k, h) - hydromet(i, km1, h)))) >
                 eps)) {

              // Ensure the abs( correlation ) between w and hydromet does not
              // have a value greater than one.
              K_hm(i, k, h) =
                  ((K_hm(i, k, h)) <
                           (((Kokkos::sqrt(wp2(i, k)) * Kokkos::sqrt(hydrometp2(i, k, h))) /
                             Kokkos::abs(
                                 (gr_invrs_dzm(i, k) * (hydromet(i, k, h) - hydromet(i, km1, h))))))
                       ? (K_hm(i, k, h))
                       : (((Kokkos::sqrt(wp2(i, k)) * Kokkos::sqrt(hydrometp2(i, k, h))) /
                           Kokkos::abs(
                               (gr_invrs_dzm(i, k) * (hydromet(i, k, h) - hydromet(i, km1, h)))))));

              // | d<hm>/dz | > 0
            }

            // Native loop: k = 1, gr%nzm-1, 1
          }

          // Set K_hm at the lower and upper boundaries (not used in calculations).
          K_hm(i, 0, h)              = zero;
          K_hm(i, ((gr_nzm - 1)), h) = zero;

          // Native loop: i = 1, hydromet_dim, 1
        }
      });

  return K_hm;
}
inline void get_cloud_top_level(int nzt, int ngrdcol, Real2D rcm, Real3D hydromet, int hydromet_dim,
                                int iiri, Int1D cloud_top_level) {
  // Description:
  // Find cloud top at a given model time step.  This function finds cloud top
  // by looping downward from the top of the model and returning the index of
  // the first vertical level that has a mean cloud water mixing ratio (or a
  // mean cloud ice mixing ratio, when ice is included in the microphysics
  // scheme) that is greater than the tolerance amount.  In a scenario that
  // there is not any cloud found, the function returns coordinate 0 (the
  // first vertical level, which is below the model surface).
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable

  //---------------------------- Begin Code ----------------------------

  // Start at the model upper boundary and loop downwards until cloud top is
  // found or the model lower boundary is reached.
  // Keep the column index innermost so each level is traversed contiguously.
  // Include cloud ice when the microphysics scheme includes ice mixing ratio.
  if ((iiri > -1)) {
    Kokkos::parallel_for(
        "get_cloud_top_level::loop2", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i        = column;
          cloud_top_level(i) = 0;
          for (k = ((nzt - 1)); k >= 0; --k) {
            if ((cloud_top_level(i) == 0)) {
              if (((rcm(i, k) > rc_tol) || (hydromet(i, k, iiri) > ri_tol))) {
                // A level with mean cloud water mixing ratio (or mean cloud ice mixing
                // ratio, when it is included in the microphysics scheme) that is
                // greater than the tolerance amount has been found.  Cloud top has
                // been found.
                cloud_top_level(i) = k;
              }
            }
          }
        });
  } else {
    Kokkos::parallel_for(
        "get_cloud_top_level::loop3", ngrdcol, KOKKOS_LAMBDA(int column) {
          int k; // Loop indices
          const int i        = column;
          cloud_top_level(i) = 0;
          for (k = ((nzt - 1)); k >= 0; --k) {
            if (((cloud_top_level(i) == 0) && (rcm(i, k) > rc_tol))) {
              cloud_top_level(i) = k;
            }
          }
        });
    // iiri > 0
  }
}
inline Int1D get_cloud_top_level(int nzt, int ngrdcol, Real2D rcm, Real3D hydromet,
                                 int hydromet_dim, int iiri) {
  // Retain owning results for callers that store cloud-top metadata persistently.
  Int1D result(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "get_cloud_top_level::cloud_top_level"),
      ngrdcol);
  get_cloud_top_level(nzt, ngrdcol, rcm, hydromet, hydromet_dim, iiri, result);
  return result;
}

inline void write_adv_micro_errors(const grid& gr, int ngrdcol, real_t dt, double time_current,
                                   int hydromet_dim, Real2D wm_zt, Real2D wp2, Real2D exner,
                                   Real2D rho, Real2D rho_zm, Real2D rcm, Real2D cloud_frac,
                                   Real2D Kh_zm, Real2D Skw_zm, Real2D rho_ds_zm, Real2D rho_ds_zt,
                                   Real2D invrs_rho_ds_zt, Real3D hydromet_mc, Real2D Ncm_mc,
                                   Real2D Lscale, Real3D hydromet_vel_covar_zt_impc,
                                   Real3D hydromet_vel_covar_zt_expc, Real2D clubb_params,
                                   const nu_vertical_res_dep& nu_vert_res_dep, bool l_upwind_xm_ma,
                                   Real3D hydromet, Real3D hydromet_vel_zt, Real3D hydrometp2,
                                   Real3D K_hm, Real2D Ncm, Real2D Nc_in_cloud, Real2D rvm_mc,
                                   Real2D thlm_mc, Real3D wphydrometp, Real2D wpNcp,
                                   const err_info_type& err_info) {
  // Description:
  // Writes to screen the values of all variables that are passed into and out
  // of subroutine advance_microphys if a fatal error has been detected.
  //-----------------------------------------------------------------------
  // Input Variables
  // approximation rather than a centered differencing for turbulent or
  // mean advection terms. It affects rtm, thlm, sclrm, um and vm.
  // Input/Output Variables for advance_microphys
  // Output Variables for advance_microphys

  //---------------------------- Begin Code ----------------------------

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {

      std::cerr << std::string("Error in advance_microphys") << "\n";

      std::cerr << std::string("Intent(in)") << "\n";

      std::cerr << std::string("dt = ") << " " << dt << "\n";
      std::cerr << std::string("time_current = ") << " " << time_current << "\n";

      std::cerr << std::string("wm_zt = ") << " " << microphys_detail::dump_values(wm_zt) << "\n";
      std::cerr << std::string("wp2 = ") << " " << microphys_detail::dump_values(wp2) << "\n";
      std::cerr << std::string("exner = ") << " " << microphys_detail::dump_values(exner) << "\n";
      std::cerr << std::string("rho = ") << " " << microphys_detail::dump_values(rho) << "\n";
      std::cerr << std::string("rho_zm = ") << " " << microphys_detail::dump_values(rho_zm) << "\n";
      std::cerr << std::string("rcm = ") << " " << microphys_detail::dump_values(rcm) << "\n";
      std::cerr << std::string("cloud_frac = ") << " " << microphys_detail::dump_values(cloud_frac)
                << "\n";
      std::cerr << std::string("Kh_zm = ") << " " << microphys_detail::dump_values(Kh_zm) << "\n";
      std::cerr << std::string("Skw_zm = ") << " " << microphys_detail::dump_values(Skw_zm) << "\n";
      std::cerr << std::string("rho_ds_zm = ") << " " << microphys_detail::dump_values(rho_ds_zm)
                << "\n";
      std::cerr << std::string("rho_ds_zt = ") << " " << microphys_detail::dump_values(rho_ds_zt)
                << "\n";
      std::cerr << std::string("invrs_rho_ds_zt = ") << " "
                << microphys_detail::dump_values(invrs_rho_ds_zt) << "\n";

      std::cerr << std::string("hydromet_mc = ") << " "
                << microphys_detail::dump_values(hydromet_mc) << "\n";

      std::cerr << std::string("Ncm_mc = ") << " " << microphys_detail::dump_values(Ncm_mc) << "\n";
      std::cerr << std::string("Lscale = ") << " " << microphys_detail::dump_values(Lscale) << "\n";

      std::cerr << std::string("hydromet_vel_covar_zt_impc = ") << " "
                << microphys_detail::dump_values(hydromet_vel_covar_zt_impc) << "\n";
      std::cerr << std::string("hydromet_vel_covar_zt_expc = ") << " "
                << microphys_detail::dump_values(hydromet_vel_covar_zt_expc) << "\n";

      std::cerr << std::string("clubb_params = ") << " "
                << microphys_detail::dump_values(clubb_params) << "\n";
      std::cerr << std::string("nu_vert_res_dep%nu_hm = ") << " "
                << microphys_detail::dump_values(nu_vert_res_dep.nu_hm) << "\n";
      std::cerr << std::string("l_upwind_xm_ma = ") << " " << l_upwind_xm_ma << "\n";

      std::cerr << std::string("Intent(inout)") << "\n";

      std::cerr << std::string("hydromet = ") << " " << microphys_detail::dump_values(hydromet)
                << "\n";
      std::cerr << std::string("hydromet_vel_zt = ") << " "
                << microphys_detail::dump_values(hydromet_vel_zt) << "\n";
      std::cerr << std::string("hydrometp2 = ") << " " << microphys_detail::dump_values(hydrometp2)
                << "\n";
      std::cerr << std::string("K_hm = ") << " " << microphys_detail::dump_values(K_hm) << "\n";

      std::cerr << std::string("Ncm = ") << " " << microphys_detail::dump_values(Ncm) << "\n";
      std::cerr << std::string("Nc_in_cloud = ") << " "
                << microphys_detail::dump_values(Nc_in_cloud) << "\n";

      std::cerr << std::string("rvm_mc = ") << " " << microphys_detail::dump_values(rvm_mc) << "\n";
      std::cerr << std::string("thlm_mc = ") << " " << microphys_detail::dump_values(thlm_mc)
                << "\n";

      std::cerr << std::string("Intent(out)") << "\n";

      std::cerr << std::string("wphydrometp = ") << " "
                << microphys_detail::dump_values(wphydrometp) << "\n";

      std::cerr << std::string("wpNcp = ") << " " << microphys_detail::dump_values(wpNcp) << "\n";
    }
  }

  return;
}
// Host/strided callers stage through the existing array owner; device Views alias.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27>
inline void
advance_microphys(ScratchWorkspace& workspace, const grid& gr, int ngrdcol, real_t dt,
                  double time_current, int hydromet_dim, const hm_metadata_type& hm_metadata,
                  Array0 wm_zt, Array1 wp2, Array2 exner, Array3 rho, Array4 rho_zm, Array5 rcm,
                  Array6 cloud_frac, Array7 Kh_zm, Array8 Skw_zm, Array9 rho_ds_zm,
                  Array10 rho_ds_zt, Array11 invrs_rho_ds_zt, Array12 hydromet_mc, Array13 Ncm_mc,
                  Array14 Lscale, Array15 hydromet_vel_covar_zt_impc,
                  Array16 hydromet_vel_covar_zt_expc, Array17 clubb_params,
                  const nu_vertical_res_dep& nu_vert_res_dep, int tridiag_solve_method,
                  int fill_holes_type, bool l_upwind_xm_ma, stats_type& stats, Array18 hydromet,
                  Array19 hydromet_vel_zt, Array20 hydrometp2, Array21 K_hm, Array22 Ncm,
                  Array23 Nc_in_cloud, Array24 rvm_mc, Array25 thlm_mc, err_info_type& err_info,
                  Array26 wphydrometp, Array27 wpNcp, const parameters_microphys& microphys) {
  // Native argument order is unchanged; only memory-space staging occurs here.
  auto workspace_scope           = workspace.scope();
  auto wm_zt_execution           = to_execution_space(workspace, wm_zt, true);
  auto wp2_execution             = to_execution_space(workspace, wp2, true);
  auto exner_execution           = to_execution_space(workspace, exner, true);
  auto rho_execution             = to_execution_space(workspace, rho, true);
  auto rho_zm_execution          = to_execution_space(workspace, rho_zm, true);
  auto rcm_execution             = to_execution_space(workspace, rcm, true);
  auto cloud_frac_execution      = to_execution_space(workspace, cloud_frac, true);
  auto Kh_zm_execution           = to_execution_space(workspace, Kh_zm, true);
  auto Skw_zm_execution          = to_execution_space(workspace, Skw_zm, true);
  auto rho_ds_zm_execution       = to_execution_space(workspace, rho_ds_zm, true);
  auto rho_ds_zt_execution       = to_execution_space(workspace, rho_ds_zt, true);
  auto invrs_rho_ds_zt_execution = to_execution_space(workspace, invrs_rho_ds_zt, true);
  auto hydromet_mc_execution     = to_execution_space(workspace, hydromet_mc, true);
  auto Ncm_mc_execution          = to_execution_space(workspace, Ncm_mc, true);
  auto Lscale_execution          = to_execution_space(workspace, Lscale, true);
  auto hydromet_vel_covar_zt_impc_execution =
      to_execution_space(workspace, hydromet_vel_covar_zt_impc, true);
  auto hydromet_vel_covar_zt_expc_execution =
      to_execution_space(workspace, hydromet_vel_covar_zt_expc, true);
  auto clubb_params_execution    = to_execution_space(workspace, clubb_params, true);
  auto hydromet_execution        = to_execution_space(workspace, hydromet, true);
  auto hydromet_vel_zt_execution = to_execution_space(workspace, hydromet_vel_zt, true);
  auto hydrometp2_execution      = to_execution_space(workspace, hydrometp2, true);
  auto K_hm_execution            = to_execution_space(workspace, K_hm, true);
  auto Ncm_execution             = to_execution_space(workspace, Ncm, true);
  auto Nc_in_cloud_execution     = to_execution_space(workspace, Nc_in_cloud, true);
  auto rvm_mc_execution          = to_execution_space(workspace, rvm_mc, true);
  auto thlm_mc_execution         = to_execution_space(workspace, thlm_mc, true);
  auto wphydrometp_execution     = to_execution_space(workspace, wphydrometp, false);
  auto wpNcp_execution           = to_execution_space(workspace, wpNcp, false);
  advance_microphys(
      workspace, gr, ngrdcol, dt, time_current, hydromet_dim, hm_metadata, wm_zt_execution,
      wp2_execution, exner_execution, rho_execution, rho_zm_execution, rcm_execution,
      cloud_frac_execution, Kh_zm_execution, Skw_zm_execution, rho_ds_zm_execution,
      rho_ds_zt_execution, invrs_rho_ds_zt_execution, hydromet_mc_execution, Ncm_mc_execution,
      Lscale_execution, hydromet_vel_covar_zt_impc_execution, hydromet_vel_covar_zt_expc_execution,
      clubb_params_execution, nu_vert_res_dep, tridiag_solve_method, fill_holes_type,
      l_upwind_xm_ma, stats, hydromet_execution, hydromet_vel_zt_execution, hydrometp2_execution,
      K_hm_execution, Ncm_execution, Nc_in_cloud_execution, rvm_mc_execution, thlm_mc_execution,
      err_info, wphydrometp_execution, wpNcp_execution, microphys);
  from_execution_space(workspace, hydromet, hydromet_execution);
  from_execution_space(workspace, hydromet_vel_zt, hydromet_vel_zt_execution);
  from_execution_space(workspace, hydrometp2, hydrometp2_execution);
  from_execution_space(workspace, K_hm, K_hm_execution);
  from_execution_space(workspace, Ncm, Ncm_execution);
  from_execution_space(workspace, Nc_in_cloud, Nc_in_cloud_execution);
  from_execution_space(workspace, rvm_mc, rvm_mc_execution);
  from_execution_space(workspace, thlm_mc, thlm_mc_execution);
  from_execution_space(workspace, wphydrometp, wphydrometp_execution);
  from_execution_space(workspace, wpNcp, wpNcp_execution);
}
} // namespace clubb_kokkos
