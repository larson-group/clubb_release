// Prepare hydrometeor PDF fields and optional SILHS sample points for microphysics.
//
// Native setup, timers, error returns and sampling stages stay in source order.
// PDF preparation, mixed moments and SILHS sampling use the C++ package APIs.
// Arrays belong to the caller. Appended configuration, permutation and RNG
// records replace implicit Native package state. The public APIs handle their
// host/device boundaries; this coordinator preserves the Native stage order.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/code_timer_module.hpp"
#include "CLUBB_core/setup_clubb_pdf_params.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "SILHS/parameters_silhs.hpp"
#include "Microphys/parameters_microphys.hpp"
#include "Microphys/mixed_moment_PDF_integrals.hpp"
#include "Microphys/microphys_diagnostics.hpp"
#include "SILHS/silhs_api_module.hpp"
#include "SILHS/latin_hypercube_driver_module.hpp"
#include "SILHS/latin_hypercube_arrays.hpp"
#include "SILHS/mt95.hpp"
#include <bit>
#include <cstdint>
#include <iostream>

namespace clubb_kokkos {
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Field7, class Field8, class Field9, class Field10, class Field11,
          class Field12, class Field13, class Field14, class Field15, class Field16, class Field17,
          class Field18, class Field19, class Field20, class Field21, class Field22, class Field23,
          class Field24, class Field25, class Field26, class Field27, class Field28, class Field29,
          class Field30>
inline void pdf_hydromet_microphys_prep(
    ScratchWorkspace& workspace,                                                        // InOut
    const grid& gr, int ngrdcol, int pdf_dim, int hydromet_dim,                         // In
    int itime, real_t vert_decorr_coef,                                                 // In
    Field0 Nc_in_cloud, Field1 cloud_frac, Field2 ice_supersat_frac,                    // In
    Field3 rho_ds_zt, Field4 Lscale, Field5 Kh_zm, Field6 hydromet, Field7 wphydrometp, // In
    Field8 corr_array_n_cloud, Field9 corr_array_n_below,                               // In
    const hm_metadata_type& hm_metadata, const pdf_parameter& pdf_params,
    Field10 clubb_params, // In
    const clubb_config_flags_type& clubb_config_flags,
    const silhs_config_flags_type& silhs_config_flags,                    // In
    bool l_rad_itime, stats_type& stats,                                  // In / InOut
    err_info_type& err_info,                                              // InOut
    Field11 hydrometp2,                                                   // Out
    Field12 mu_x_1_n, Field13 mu_x_2_n,                                   // Out
    Field14 sigma_x_1_n, Field15 sigma_x_2_n,                             // Out
    Field16 corr_array_1_n, Field17 corr_array_2_n,                       // Out
    Field18 corr_cholesky_mtx_1, Field19 corr_cholesky_mtx_2,             // Out
    precipitation_fractions& precip_fracs,                                // InOut
    Field20 rtphmp_zt, Field21 thlphmp_zt, Field22 wp2hmp,                // Out
    Field23 X_nl_all_levs, Field24 X_mixt_comp_all_levs,                  // Out
    Field25 lh_sample_point_weights,                                      // Out
    Field26 lh_rt_clipped, Field27 lh_thl_clipped, Field28 lh_rc_clipped, // Out
    Field29 lh_rv_clipped, Field30 lh_Nc_clipped,                         // Out
    hydromet_pdf_parameter_collection* hydromet_pdf_params,               // Optional Out, nullable
    const parameters_microphys& microphys,        // Case-owned Native microphysics settings.
    const parameters_silhs_state& silhs,          // Case-owned Native importance settings.
    latin_hypercube_arrays_state& sampling_state, // Case-owned permutation/cache state.
    mt95& rng                                     // Case-owned Native Mersenne Twister state.
) {
  // Description:
  // Prepares the hydrometeor PDF and, when SILHS is in use, sample points,
  // for the call to microphysics. Outputs retain the Native inactive/early-
  // return definedness; no JAX-only initialization is added.
  // gr: Grid metadata and coordinate fields.
  // ngrdcol: Number of grid columns
  // pdf_dim: Number of variables in the correlation array
  // hydromet_dim: Number of hydrometeor species
  // itime: Current model timestep index.
  // vert_decorr_coef: Empirically defined de-correlation constant [-]
  // Nc_in_cloud: Mean (in-cloud) cloud droplet conc.       [num/kg]
  // cloud_frac: Cloud fraction                            [-]
  // ice_supersat_frac: Ice supersaturation fraction              [-]
  // rho_ds_zt: Dry, base-state density on thermo. levs.  [kg/m^3]
  // Lscale: Turbulent Mixing Length                   [m]
  // Kh_zm: Eddy diffusivity coef. on momentum levels [m^2/s]
  // hydromet: Mean of hydrometeor, hm (overall) (t-levs.) [units]
  // wphydrometp: Covariance < w'h_m' > (momentum levels)     [(m/s)units]
  // corr_array_n_cloud: Prescribed normal space corr. array in cloud  [-]
  // corr_array_n_below: Prescribed normal space corr. array below cl. [-]
  // hm_metadata: Hydrometeor/PDF indices, names, flags and tolerances.
  // pdf_params: PDF parameters                               [units vary]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // clubb_config_flags: Derived type holding all configurable CLUBB flags
  // silhs_config_flags: SILHS sampling configuration.
  // l_rad_itime: Whether radiation is called at this timestep.
  // stats: Native statistics collector; may accumulate output.
  // err_info: err_info struct containing err_code and err_header
  // hydrometp2: Variance of a hydrometeor (overall) (m-levs.)   [units^2]
  // mu_x_1_n: Mean array (normal space): PDF vars. (comp. 1) [un. vary]
  // mu_x_2_n: Mean array (normal space): PDF vars. (comp. 2) [un. vary]
  // sigma_x_1_n: Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
  // sigma_x_2_n: Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]
  // corr_array_1_n: Corr. array (normal space) of PDF vars. (comp. 1)  [-]
  // corr_array_2_n: Corr. array (normal space) of PDF vars. (comp. 2)  [-]
  // corr_cholesky_mtx_1: Transposed corr. cholesky matrix, 1st comp. [-]
  // corr_cholesky_mtx_2: Transposed corr. cholesky matrix, 2nd comp. [-]
  // precip_fracs: Precipitation fractions      [-]
  // rtphmp_zt: Covariance of rt and hm (on t-levs.)     [(kg/kg)<hm un.>]
  // thlphmp_zt: Covariance of thl and hm (on t-levs.)    [K<hm units>]
  // wp2hmp: Higher-order mixed moment:  < w'^2 hm' > [(m/s)^2<hm un.>]
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // X_mixt_comp_all_levs: Which mixture component we're in
  // lh_sample_point_weights: Sample-point weights for SILHS statistics [-].
  // lh_rt_clipped: rt generated from silhs sample points
  // lh_thl_clipped: thl generated from silhs sample points
  // lh_rc_clipped: rc generated from silhs sample points
  // lh_rv_clipped: rv generated from silhs sample points
  // lh_Nc_clipped: Nc generated from silhs sample points
  // hydromet_pdf_params: Hydrometeor PDF parameters        [units vary]
  // microphys: Native microphysics settings used for scheme/count/seed selection.
  // silhs: Native prescribed category probabilities and importance threshold.
  // sampling_state: Saved Native one_height_time_matrix and prior iteration.
  // rng: Native generator state; sampler resets/advances it in source order.

  // Local Variables
  auto workspace_scope = workspace.scope();
  bool l_calc_weights_all_levs_itime; // Whether every level needs distinct sample weights.
  constexpr bool l_calc_weights_all_levs = false; // Same weights at every grid level.
  // Whether to call Ncn_to_Nc. It may cause MG problems because the resulting
  // Nc tendency is not fed back into the microphysics, as noted in Native.
  constexpr bool l_use_Ncn_to_Nc = true;
  std::int32_t lh_seed_custom; // Native genrand_intg=selected_int_kind(9), the current 32-bit kind.

  // Required state follows the full Native prefix. The optional output remains
  // nullable; callers pass nullptr explicitly because required state follows it.

  // -------------------------------- Begin Code --------------------------------
  if (microphys_detail::trim(microphys.microphys_scheme) != "none") {
    // Host/device updates belong to the called APIs and their existing array owner.
    // Setup the PDF parameters.
    timer_start("setup_pdf_parameters_api");
    setup_pdf_parameters_api(workspace, gr, gr.nzm, gr.nzt, ngrdcol, pdf_dim, hydromet_dim, // In
                             Nc_in_cloud, cloud_frac, Kh_zm,                                // In
                             ice_supersat_frac, hydromet, wphydrometp,                      // In
                             corr_array_n_cloud, corr_array_n_below,                        // In
                             hm_metadata,                                                   // In
                             pdf_params,                                                    // In
                             clubb_params,                                                  // In
                             clubb_config_flags.iiPDF_type,                                 // In
                             clubb_config_flags.l_use_precip_frac,                          // In
                             clubb_config_flags.l_diagnose_correlations,                    // In
                             clubb_config_flags.l_calc_w_corr,                              // In
                             clubb_config_flags.l_const_Nc_in_cloud,                        // In
                             clubb_config_flags.l_fix_w_chi_eta_correlations,               // In
                             err_info,                                                      // InOut
                             hydrometp2,                                                    // Out
                             mu_x_1_n, mu_x_2_n,                                            // Out
                             sigma_x_1_n, sigma_x_2_n,                                      // Out
                             corr_array_1_n, corr_array_2_n,                                // Out
                             corr_cholesky_mtx_1, corr_cholesky_mtx_2,                      // Out
                             precip_fracs,                                                  // InOut
                             hydromet_pdf_params, // Optional Out
                             stats);              // InOut
    timer_stop("setup_pdf_parameters_api");

    // Error check after setup_pdf_parameters.
    if (clubb_at_least_debug_level_api(0)) {
      auto codes = to_host_space(workspace, err_info.err_code);
      bool fatal = false;
      for (std::size_t i = 0; i < codes.extent(0); ++i)
        fatal = fatal || codes(i) == clubb_fatal_error;
      if (fatal) {
        std::cerr << err_info.err_header_global << '\n';
        std::cerr << "Fatal error after setup_pdf_parameters_api" << '\n';
        return;
      }
    }

    // Calculate <rt'hm'>, <thl'hm'> and <w'^2 hm'>.
    timer_start("hydrometeor_mixed_moments");
    // Native forwards the optional collection to this required lower input.
    // An active microphysics call therefore requires its presence; the driver
    // supplies it. Absent-active behavior has no defined Native reference.
    clubb_kokkos::hydrometeor_mixed_moments(workspace, gr, ngrdcol, gr.nzt, pdf_dim,
                                            hydromet_dim,                     // In
                                            hydromet, hm_metadata,            // In
                                            mu_x_1_n, mu_x_2_n,               // In
                                            sigma_x_1_n, sigma_x_2_n,         // In
                                            corr_array_1_n, corr_array_2_n,   // In
                                            pdf_params, *hydromet_pdf_params, // In
                                            precip_fracs, stats,              // In / InOut
                                            rtphmp_zt, thlphmp_zt, wp2hmp);   // Out
    timer_stop("hydrometeor_mixed_moments");
  }

  //----------------------------------------------------------------
  // Compute subcolumns if enabled.
  //----------------------------------------------------------------
  if (microphys.lh_microphys_type != lh_microphys_disabled || l_silhs_rad) {
    // Calculate sample weights separately at all grid levels when radiation is not called.
    l_calc_weights_all_levs_itime = l_calc_weights_all_levs && !l_rad_itime;

    // The profile of CLUBB's mixing length, Lscale, is passed into
    // subroutine generate_silhs_sample for use in calculating the vertical
    // correlation coefficient.  Rather than output Lscale directly, its
    // value can be calculated from other fields that are already output from
    // subroutine advance_clubb_core.  The equation relating Lscale to eddy
    // diffusivity is:
    //
    // Kh = c_K * Lscale * sqrt( TKE ).
    //
    // Kh is available, TKE can be calculated from wp2, up2, and vp2, all of
    // which are available, and c_K is easily extracted from CLUBB's tunable
    // parameters.  The equation for Lscale is:
    //
    // Lscale = Kh / ( c_K * sqrt( TKE ) ).
    //
    // Since Kh and TKE are output on momentum grid levels, the resulting
    // calculation of Lscale is also found on momentum levels.  It needs to
    // be interpolated back to thermodynamic (midpoint) grid levels for
    // further use.

    // Rather than using the seed just to intialize the random number generator at the first
    // timestep, we use it every
    // timestep. Thus this seed defines the randoms, that way we can reproduce results for restart
    // runs by using a
    // reproducible seed. If we didn't do it this way, we would have to save the internal state of
    // the random number
    // generator to disk each timestep, then read it back in to reinitialize the random number
    // generator
    // upon restart - which is too much work to implement
    // BFB-022: Native genrand_intg is signed32. Preserve its seed product bits
    // without invoking C++ signed-overflow undefined behavior on long runs.
    lh_seed_custom = std::bit_cast<std::int32_t>(static_cast<std::uint32_t>(microphys.lh_seed) *
                                                 static_cast<std::uint32_t>(itime));

    timer_start("generate_silhs_sample_api");
    clubb_kokkos::generate_silhs_sample_api(workspace, itime, pdf_dim, microphys.lh_num_samples,
                                            microphys.lh_sequence_length, gr.nzt, ngrdcol, // In
                                            l_calc_weights_all_levs_itime,                 // In
                                            gr, pdf_params, gr.dzt, Lscale,                // In
                                            lh_seed_custom, hm_metadata,                   // In
                                            mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,  // In
                                            corr_cholesky_mtx_1, corr_cholesky_mtx_2,      // In
                                            precip_fracs, silhs_config_flags,              // In
                                            vert_decorr_coef,                              // In
                                            err_info,                                      // InOut
                                            X_nl_all_levs, X_mixt_comp_all_levs,           // Out
                                            lh_sample_point_weights,                       // Out
                                            &stats, // Optional InOut, present
                                            sampling_state, rng,
                                            silhs); // Appended case-owned Native state
    timer_stop("generate_silhs_sample_api");

    // Error check after generating sample points.
    if (clubb_at_least_debug_level_api(0)) {
      auto codes = to_host_space(workspace, err_info.err_code);
      bool fatal = false;
      for (std::size_t i = 0; i < codes.extent(0); ++i)
        fatal = fatal || codes(i) == clubb_fatal_error;
      if (fatal) {
        std::cerr << err_info.err_header_global << '\n';
        std::cerr << "Fatal error calling generate_silhs_sample_api in "
                     "pdf_hydromet_microphys_prep"
                  << '\n';
        return;
      }
    }

    timer_start("clip_transform_silhs_output_api");
    clubb_kokkos::clip_transform_silhs_output_api(workspace, gr.nzt, ngrdcol,
                                                  microphys.lh_num_samples,           // In
                                                  pdf_dim, hydromet_dim, hm_metadata, // In
                                                  X_mixt_comp_all_levs,               // In
                                                  X_nl_all_levs,                      // InOut
                                                  pdf_params, l_use_Ncn_to_Nc,        // In
                                                  lh_rt_clipped, lh_thl_clipped,      // Out
                                                  lh_rc_clipped, lh_rv_clipped,       // Out
                                                  lh_Nc_clipped);                     // Out
    timer_stop("clip_transform_silhs_output_api");

    if (stats.l_sample) {
      timer_start("stats_accumulate_lh_api");
      clubb_kokkos::stats_accumulate_lh_api(workspace, gr, gr.nzt, ngrdcol,
                                            microphys.lh_num_samples, pdf_dim,
                                            rho_ds_zt,                                   // In
                                            hydromet_dim, hm_metadata,                   // In
                                            lh_sample_point_weights, X_nl_all_levs,      // In
                                            lh_rt_clipped, lh_thl_clipped,               // In
                                            lh_rc_clipped, lh_rv_clipped, lh_Nc_clipped, // In
                                            &stats); // Optional InOut, present
      timer_stop("stats_accumulate_lh_api");
    }
  }
}
} // namespace clubb_kokkos
