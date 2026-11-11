// Generate, transform and diagnose SILHS sample points in all columns.
// Public functions: generate_silhs_sample, clip_transform_silhs_output,
// latin_hypercube_2D_output_api, stats_accumulate_lh_api and
// copy_X_nl_into_hydromet_all_pts. State and optional C++ stats are case owned.
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "SILHS/parameters_silhs.hpp"
#include "SILHS/latin_hypercube_arrays.hpp"
#include "SILHS/mt95.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

namespace clubb_kokkos {
void generate_silhs_sample(ScratchWorkspace& workspace, int iter, int pdf_dim, int num_samples,
                           int sequence_length, int nzt,
                           int ngrdcol,                        // In
                           bool l_calc_weights_all_levs_itime, // In
                           const grid& gr, const pdf_parameter& pdf_params, Real2D delta_zm,
                           Real2D Lscale,                                          // In
                           int lh_seed, const hm_metadata_type& hm_metadata,       // In
                           Real3D mu1, Real3D mu2, Real3D sigma1, Real3D sigma2,   // In
                           Real4D corr_cholesky_mtx_1, Real4D corr_cholesky_mtx_2, // In
                           const precipitation_fractions& precip_fracs,            // In
                           const silhs_config_flags_type& silhs_config_flags,      // In
                           real_t vert_decorr_coef,                                // In
                           err_info_type& err_info,                                // InOut
                           Real4D X_nl_all_levs,
                           Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs, // Out
                           Real3D lh_sample_point_weights,                                // Out
                           stats_type* stats,                            // Optional InOut
                           latin_hypercube_arrays_state& sampling_state, // Case state
                           mt95& rng,                                    // Case state
                           const parameters_silhs_state& silhs           // Case state
);

void clip_transform_silhs_output(
    ScratchWorkspace& workspace,                                        // InOut
    int nzt, int ngrdcol, int num_samples,                              // In
    int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata, // In
    Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,      // In
    Real4D X_nl_all_levs,                                               // InOut
    const pdf_parameter& pdf_params, bool l_use_Ncn_to_Nc,              // In
    Real3D lh_rt_clipped, Real3D lh_thl_clipped,                        // Out
    Real3D lh_rc_clipped, Real3D lh_rv_clipped,                         // Out
    Real3D lh_Nc_clipped                                                // Out
);

void latin_hypercube_2D_output_api(int nzt,
                                   Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> zt,
                                   int pdf_dim, int num_samples,        // In
                                   const hm_metadata_type& hm_metadata, // In
                                   stats_type* stats,                   // Optional InOut
                                   err_info_type& err_info              // InOut
);

void stats_accumulate_lh_api(ScratchWorkspace& workspace, const grid& gr, int nzt, int ngrdcol,
                             int num_samples, int pdf_dim,
                             Real2D rho_ds_zt,                                      // In
                             int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                             Real3D lh_sample_point_weights, Real4D X_nl_all_levs,  // In
                             Real3D lh_rt_clipped, Real3D lh_thl_clipped,           // In
                             Real3D lh_rc_clipped, Real3D lh_rv_clipped,            // In
                             Real3D lh_Nc_clipped,                                  // In
                             stats_type* stats // Optional InOut
);

void copy_X_nl_into_hydromet_all_pts(int nzt, int pdf_dim, int num_samples, int ngrdcol,
                                     Real4D X_nl_all_levs, // In
                                     int hydromet_dim, const hm_metadata_type& hm_metadata,
                                     Real3D hydromet,                                  // In
                                     Real4D hydromet_all_points, Real3D Ncn_all_points // Out
);

template <class Field_delta_zm, class Field_Lscale, class Field_mu1, class Field_mu2,
          class Field_sigma1, class Field_sigma2, class Field_corr_cholesky_mtx_1,
          class Field_corr_cholesky_mtx_2, class Field_X_nl_all_levs,
          class Field_X_mixt_comp_all_levs, class Field_lh_sample_point_weights>
void generate_silhs_sample(
    ScratchWorkspace& workspace,                                                        // InOut
    int iter, int pdf_dim, int num_samples, int sequence_length, int nzt, int ngrdcol,  // In
    bool l_calc_weights_all_levs_itime,                                                 // In
    const grid& gr, const pdf_parameter& pdf_params, Field_delta_zm delta_zm,           // In
    Field_Lscale Lscale,                                                                // In
    int lh_seed, const hm_metadata_type& hm_metadata,                                   // In
    Field_mu1 mu1, Field_mu2 mu2, Field_sigma1 sigma1, Field_sigma2 sigma2,             // In
    Field_corr_cholesky_mtx_1 corr_cholesky_mtx_1,                                      // In
    Field_corr_cholesky_mtx_2 corr_cholesky_mtx_2,                                      // In
    const precipitation_fractions& precip_fracs,                                        // In
    const silhs_config_flags_type& silhs_config_flags,                                  // In
    real_t vert_decorr_coef,                                                            // In
    err_info_type& err_info,                                                            // InOut
    Field_X_nl_all_levs X_nl_all_levs, Field_X_mixt_comp_all_levs X_mixt_comp_all_levs, // Out
    Field_lh_sample_point_weights lh_sample_point_weights,                              // Out
    stats_type* stats,                            // Optional InOut
    latin_hypercube_arrays_state& sampling_state, // Case state
    mt95& rng,                                    // Case state
    const parameters_silhs_state& silhs           // Case state
) {
  // Shared-owner staging at the public host/strided field boundary.
  // iter: Model iteration (time step) number
  // pdf_dim: Number of variables to sample
  // num_samples: Number of samples per variable
  // sequence_length: nt_repeat/num_samples; number of timesteps before sequence repeats
  // nzt: Number of thermodynamic vertical model levels
  // ngrdcol: Number of grid columns
  // l_calc_weights_all_levs_itime: determines if vertically correlated sample points are needed
  // gr: Grid variable type
  // pdf_params: PDF parameters       [units vary]
  // delta_zm: Difference in momentum altitudes    [m]
  // Lscale: Turbulent mixing length            [m]
  // lh_seed: Random number generator seed
  // mu1: Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
  // mu2: Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
  // sigma1: Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
  // sigma2: Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]
  // corr_cholesky_mtx_1: Correlations Cholesky matrix (1st comp.)  [-]
  // corr_cholesky_mtx_2: Correlations Cholesky matrix (2nd comp.)  [-]
  // precip_fracs: Precipitation fractions      [-]
  // silhs_config_flags: Flags for the SILHS sampling code [-]
  // vert_decorr_coef: Empirically defined de-correlation constant [-]
  // err_info: err_info struct containing err_code and err_header
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // X_mixt_comp_all_levs: Which mixture component we're in
  // lh_sample_point_weights: Weight of each sample point
  auto workspace_scope                = workspace.scope();
  auto X_nl_all_levs_execution        = to_execution_space(workspace, X_nl_all_levs, false);
  auto X_mixt_comp_all_levs_execution = to_execution_space(workspace, X_mixt_comp_all_levs, false);
  auto lh_sample_point_weights_execution =
      to_execution_space(workspace, lh_sample_point_weights, false);
  generate_silhs_sample(workspace, iter, pdf_dim, num_samples, sequence_length, nzt, ngrdcol, // In
                        l_calc_weights_all_levs_itime,                                        // In
                        gr, pdf_params, to_execution_space(workspace, delta_zm),
                        to_execution_space(workspace, Lscale), // In
                        lh_seed, hm_metadata,                  // In
                        to_execution_space(workspace, mu1), to_execution_space(workspace, mu2),
                        to_execution_space(workspace, sigma1),
                        to_execution_space(workspace, sigma2), // In
                        to_execution_space(workspace, corr_cholesky_mtx_1),
                        to_execution_space(workspace, corr_cholesky_mtx_2),      // In
                        precip_fracs, silhs_config_flags,                        // In
                        vert_decorr_coef,                                        // In
                        err_info,                                                // InOut
                        X_nl_all_levs_execution, X_mixt_comp_all_levs_execution, // Out
                        lh_sample_point_weights_execution,                       // Out
                        stats,                                                   // Optional InOut
                        sampling_state, rng, silhs                               // Case state
  );
  from_execution_space(workspace, X_nl_all_levs, X_nl_all_levs_execution);
  from_execution_space(workspace, X_mixt_comp_all_levs, X_mixt_comp_all_levs_execution);
  from_execution_space(workspace, lh_sample_point_weights, lh_sample_point_weights_execution);
}

template <class Field_X_mixt_comp_all_levs, class Field_X_nl_all_levs, class Field_lh_rt_clipped,
          class Field_lh_thl_clipped, class Field_lh_rc_clipped, class Field_lh_rv_clipped,
          class Field_lh_Nc_clipped>
void clip_transform_silhs_output(ScratchWorkspace& workspace,           // InOut
                                 int nzt, int ngrdcol, int num_samples, // In
                                 int pdf_dim, int hydromet_dim,
                                 const hm_metadata_type& hm_metadata,                   // In
                                 Field_X_mixt_comp_all_levs X_mixt_comp_all_levs,       // In
                                 Field_X_nl_all_levs X_nl_all_levs,                     // InOut
                                 const pdf_parameter& pdf_params, bool l_use_Ncn_to_Nc, // In
                                 Field_lh_rt_clipped lh_rt_clipped,
                                 Field_lh_thl_clipped lh_thl_clipped, // Out
                                 Field_lh_rc_clipped lh_rc_clipped,
                                 Field_lh_rv_clipped lh_rv_clipped, // Out
                                 Field_lh_Nc_clipped lh_Nc_clipped  // Out
) {
  // Shared-owner staging at the public host/strided field boundary.
  // nzt: Number of vertical levels
  // ngrdcol: Number of grid columns
  // num_samples: Number of SILHS sample points
  // pdf_dim: Number of variates in X_nl
  // hydromet_dim: Number of hydrometeor species
  // X_mixt_comp_all_levs: Which component this sample is in (1 or 2)
  // X_nl_all_levs: SILHS sample points    [units vary]
  // pdf_params: **The** PDF parameters!
  // l_use_Ncn_to_Nc: Whether to call Ncn_to_Nc (.true.) or not (.false.); Ncn_to_Nc might cause
  //   problems with the MG microphysics since the changes made here (Nc-tendency) are not fed into
  //   the microphysics
  // lh_rt_clipped: rt generated from silhs sample points
  // lh_thl_clipped: thl generated from silhs sample points
  // lh_rc_clipped: rc generated from silhs sample points
  // lh_rv_clipped: rv generated from silhs sample points
  // lh_Nc_clipped: Nc generated from silhs sample points
  auto workspace_scope          = workspace.scope();
  auto X_nl_all_levs_execution  = to_execution_space(workspace, X_nl_all_levs);
  auto lh_rt_clipped_execution  = to_execution_space(workspace, lh_rt_clipped, false);
  auto lh_thl_clipped_execution = to_execution_space(workspace, lh_thl_clipped, false);
  auto lh_rc_clipped_execution  = to_execution_space(workspace, lh_rc_clipped, false);
  auto lh_rv_clipped_execution  = to_execution_space(workspace, lh_rv_clipped, false);
  auto lh_Nc_clipped_execution  = to_execution_space(workspace, lh_Nc_clipped, false);
  clip_transform_silhs_output(workspace, nzt, ngrdcol, num_samples,                // In
                              pdf_dim, hydromet_dim, hm_metadata,                  // In
                              to_execution_space(workspace, X_mixt_comp_all_levs), // In
                              X_nl_all_levs_execution,                             // InOut
                              pdf_params, l_use_Ncn_to_Nc,                         // In
                              lh_rt_clipped_execution, lh_thl_clipped_execution,   // Out
                              lh_rc_clipped_execution, lh_rv_clipped_execution,    // Out
                              lh_Nc_clipped_execution                              // Out
  );
  from_execution_space(workspace, X_nl_all_levs, X_nl_all_levs_execution);
  from_execution_space(workspace, lh_rt_clipped, lh_rt_clipped_execution);
  from_execution_space(workspace, lh_thl_clipped, lh_thl_clipped_execution);
  from_execution_space(workspace, lh_rc_clipped, lh_rc_clipped_execution);
  from_execution_space(workspace, lh_rv_clipped, lh_rv_clipped_execution);
  from_execution_space(workspace, lh_Nc_clipped, lh_Nc_clipped_execution);
}

template <class Field_zt>
void latin_hypercube_2D_output_api(ScratchWorkspace& workspace, int nzt, Field_zt zt, int pdf_dim,
                                   int num_samples,                     // In
                                   const hm_metadata_type& hm_metadata, // In
                                   stats_type* stats,                   // Optional InOut
                                   err_info_type& err_info              // InOut
) {
  // SILHS file setup consumes an independent host height snapshot.
  // Return before staging when the optional stats context is disabled.
  if (!stats || !stats->enabled)
    return;
  auto workspace_scope = workspace.scope();
  latin_hypercube_2D_output_api(nzt, to_host_space(workspace, zt), pdf_dim, num_samples,
                                hm_metadata,    // In
                                stats, err_info // Optional InOut / InOut
  );
}

template <class Field_rho_ds_zt, class Field_lh_sample_point_weights, class Field_X_nl_all_levs,
          class Field_lh_rt_clipped, class Field_lh_thl_clipped, class Field_lh_rc_clipped,
          class Field_lh_rv_clipped, class Field_lh_Nc_clipped>
void stats_accumulate_lh_api(ScratchWorkspace& workspace, const grid& gr, int nzt, int ngrdcol,
                             int num_samples,
                             int pdf_dim,                                           // In
                             Field_rho_ds_zt rho_ds_zt,                             // In
                             int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                             Field_lh_sample_point_weights lh_sample_point_weights, // In
                             Field_X_nl_all_levs X_nl_all_levs,                     // In
                             Field_lh_rt_clipped lh_rt_clipped,
                             Field_lh_thl_clipped lh_thl_clipped, // In
                             Field_lh_rc_clipped lh_rc_clipped,
                             Field_lh_rv_clipped lh_rv_clipped, // In
                             Field_lh_Nc_clipped lh_Nc_clipped, // In
                             stats_type* stats                  // Optional InOut
) {
  // Shared-owner staging at the public host/strided field boundary.
  // nzt: Number of vertical model levels
  // ngrdcol: Number of model columns
  // num_samples: Number of calls to microphysics per timestep (normally=2)
  // pdf_dim: Number of variables to sample
  // rho_ds_zt: Dry, static density (thermo. levs.) [kg/m^3]
  // hydromet_dim: Number of hydrometeor species
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // lh_rt_clipped: rt generated from silhs sample points
  // lh_thl_clipped: thl generated from silhs sample points
  // lh_rc_clipped: rc generated from silhs sample points
  // lh_rv_clipped: rv generated from silhs sample points
  // lh_Nc_clipped: Nc generated from silhs sample points
  auto workspace_scope = workspace.scope();
  stats_accumulate_lh_api(workspace, gr, nzt, ngrdcol, num_samples, pdf_dim,
                          to_execution_space(workspace, rho_ds_zt), // In
                          hydromet_dim, hm_metadata,                // In
                          to_execution_space(workspace, lh_sample_point_weights),
                          to_execution_space(workspace, X_nl_all_levs), // In
                          to_execution_space(workspace, lh_rt_clipped),
                          to_execution_space(workspace, lh_thl_clipped), // In
                          to_execution_space(workspace, lh_rc_clipped),
                          to_execution_space(workspace, lh_rv_clipped), // In
                          to_execution_space(workspace, lh_Nc_clipped), // In
                          stats                                         // Optional InOut
  );
}

template <class Field_X_nl_all_levs, class Field_hydromet, class Field_hydromet_all_points,
          class Field_Ncn_all_points>
void copy_X_nl_into_hydromet_all_pts(int nzt, int pdf_dim, int num_samples, int ngrdcol, // In
                                     Field_X_nl_all_levs X_nl_all_levs,                  // In
                                     int hydromet_dim, const hm_metadata_type& hm_metadata,
                                     Field_hydromet hydromet, // In
                                     Field_hydromet_all_points hydromet_all_points,
                                     Field_Ncn_all_points Ncn_all_points // Out
) {
  // Shared-owner staging at the public host/strided field boundary.
  // hydromet_dim: Number of hydrometeor species
  // hydromet: Hydrometeor species    [units vary]
  auto hydromet_all_points_execution = to_execution_space(hydromet_all_points, false);
  auto Ncn_all_points_execution      = to_execution_space(Ncn_all_points, false);
  copy_X_nl_into_hydromet_all_pts(nzt, pdf_dim, num_samples, ngrdcol,
                                  to_execution_space(X_nl_all_levs),                       // In
                                  hydromet_dim, hm_metadata, to_execution_space(hydromet), // In
                                  hydromet_all_points_execution, Ncn_all_points_execution  // Out
  );
  from_execution_space(hydromet_all_points, hydromet_all_points_execution);
  from_execution_space(Ncn_all_points, Ncn_all_points_execution);
}

} // namespace clubb_kokkos
