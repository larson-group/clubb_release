// Public SILHS generation and clipping API. Kokkos owns array placement;
// case-owned permutation/RNG/configuration replace Native implicit state.
#include "scratch_workspace.hpp"
#include "SILHS/silhs_api_module.hpp"
#include "SILHS/latin_hypercube_driver_module.hpp"
namespace clubb_kokkos {

//-----------------------------------------------------------------------
void generate_silhs_sample_api(ScratchWorkspace& workspace, int iter, int pdf_dim, int num_samples,
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
) {
  // Generate samples in all columns; the public array adapter owns the
  // Native data-region boundary. Preserve optional statistics presence.
  // iter: Model iteration number
  // pdf_dim: Number of variables to sample
  // num_samples: Number of samples per variable
  // sequence_length: nt_repeat/num_samples; number of timesteps before sequence repeats.
  // nzt: Number of vertical model levels
  // ngrdcol: Number of grid columns
  // l_calc_weights_all_levs_itime: determines if vertically correlated sample points are needed
  // gr: Grid variable type
  // pdf_params: PDF parameters       [units vary]
  // delta_zm: Difference in moment. altitudes    [m]
  // Lscale: Turbulent Mixing Length  [m]
  // lh_seed: Random number generator seed
  // mu1: Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
  // mu2: Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
  // sigma1: Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
  // sigma2: Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]
  // corr_cholesky_mtx_1: Correlations Cholesky matrix (1st comp.)  [-]
  // corr_cholesky_mtx_2: Correlations Cholesky matrix (2nd comp.)  [-]
  // precip_fracs: Precipitation fractions      [-]
  // vert_decorr_coef: Empirically defined de-correlation constant [-]
  // err_info: err_info struct containing err_code and err_header
  // X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
  // X_mixt_comp_all_levs: Which mixture component we're in
  auto workspace_scope = workspace.scope();
  if (stats) {
    generate_silhs_sample(workspace, iter, pdf_dim, num_samples, sequence_length, nzt,
                          ngrdcol,                                  // In
                          l_calc_weights_all_levs_itime,            // In
                          gr, pdf_params, delta_zm, Lscale,         // In
                          lh_seed, hm_metadata,                     // In
                          mu1, mu2, sigma1, sigma2,                 // In
                          corr_cholesky_mtx_1, corr_cholesky_mtx_2, // In
                          precip_fracs, silhs_config_flags,         // In
                          vert_decorr_coef,                         // In
                          err_info,                                 // InOut
                          X_nl_all_levs, X_mixt_comp_all_levs,      // Out
                          lh_sample_point_weights,                  // Out
                          stats,                                    // Optional InOut
                          sampling_state, rng, silhs                // Case state
    );
  } else {
    generate_silhs_sample(workspace, iter, pdf_dim, num_samples, sequence_length, nzt,
                          ngrdcol,                                  // In
                          l_calc_weights_all_levs_itime,            // In
                          gr, pdf_params, delta_zm, Lscale,         // In
                          lh_seed, hm_metadata,                     // In
                          mu1, mu2, sigma1, sigma2,                 // In
                          corr_cholesky_mtx_1, corr_cholesky_mtx_2, // In
                          precip_fracs, silhs_config_flags,         // In
                          vert_decorr_coef,                         // In
                          err_info,                                 // InOut
                          X_nl_all_levs, X_mixt_comp_all_levs,      // Out
                          lh_sample_point_weights,                  // Out
                          nullptr,                                  // Optional InOut
                          sampling_state, rng, silhs                // Case state
    );
  }
}

//-----------------------------------------------------------------------
void clip_transform_silhs_output_api(
    ScratchWorkspace& workspace,                                        // InOut
    int nzt, int ngrdcol, int num_samples,                              // In
    int pdf_dim, int hydromet_dim, const hm_metadata_type& hm_metadata, // In
    Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,      // In
    Real4D X_nl_all_levs,                                               // InOut
    const pdf_parameter& pdf_params, bool l_use_Ncn_to_Nc,              // In
    Real3D lh_rt_clipped, Real3D lh_thl_clipped,                        // Out
    Real3D lh_rc_clipped, Real3D lh_rv_clipped,                         // Out
    Real3D lh_Nc_clipped                                                // Out
) {
  // Compute clipped rt,thl,rc,rv,Nc for every sample and column.
  // nzt: Number of vertical levels
  // ngrdcol: Number of grid columns
  // num_samples: Number of SILHS sample points
  // pdf_dim: Number of variates in X_nl_one_lev
  // hydromet_dim: Number of hydrometeor species
  // X_mixt_comp_all_levs: Which component this sample is in (1 or 2)
  // X_nl_all_levs: SILHS sample points    [units vary]
  // pdf_params: **The** PDF parameters!
  // lh_rt_clipped: rt generated from silhs sample points
  // lh_thl_clipped: thl generated from silhs sample points
  // lh_rc_clipped: rc generated from silhs sample points
  // lh_rv_clipped: rv generated from silhs sample points
  // lh_Nc_clipped: Nc generated from silhs sample points
  auto workspace_scope = workspace.scope();
  clip_transform_silhs_output(workspace, nzt, ngrdcol, num_samples, // In
                              pdf_dim, hydromet_dim, hm_metadata,   // In
                              X_mixt_comp_all_levs,                 // In
                              X_nl_all_levs,                        // InOut
                              pdf_params, l_use_Ncn_to_Nc,          // In
                              lh_rt_clipped, lh_thl_clipped,        // Out
                              lh_rc_clipped, lh_rv_clipped,         // Out
                              lh_Nc_clipped                         // Out
  );
}
} // namespace clubb_kokkos
