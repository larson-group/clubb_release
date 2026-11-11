// SILHS sampling, clipping and diagnostic statistics in Native source order.
// Physical fields use column-first execution-space Views and zero-based level
// values. Saved k_lh_start statistics retain Native ordinals. Ordered random draws
// and permutation/importance control share explicit case-owned state.
#include "scratch_workspace.hpp"
#include "SILHS/latin_hypercube_driver_module.hpp"
#include "SILHS/generate_uniform_sample_module.hpp"
#include "SILHS/silhs_importance_sample_module.hpp"
#include "SILHS/transform_to_pdf_module.hpp"
#include "SILHS/math_utilities.hpp"
#include "CLUBB_core/advance_helper_module.hpp"
#include "CLUBB_core/fill_holes.hpp"
#include "CLUBB_core/index_mapping.hpp"
#include "CLUBB_core/pdf_utilities.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/stats_netcdf.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <vector>

namespace clubb_kokkos {
namespace {
constexpr bool l_output_2D_lognormal_dist = false;
constexpr bool l_output_2D_uniform_dist   = false;
void generate_random_pool(ScratchWorkspace& workspace, int nzt, int ngrdcol, int pdf_dim,
                          int num_samples,
                          int d_uniform_extra,          // In
                          int lh_seed, const grid& gr,  // In
                          bool l_lh_deterministic_test, // In
                          Real4D rand_pool,             // Out
                          mt95& rng                     // Case state
);
void generate_all_uniform_samples(ScratchWorkspace& workspace, int iter, int pdf_dim,
                                  int d_uniform_extra, int num_samples,
                                  int sequence_length, // In
                                  int nzt, int ngrdcol,
                                  Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // In
                                  Real2D X_vert_corr, Real4D rand_pool,              // In
                                  int iiPDF_chi,                                     // In
                                  Real2D cloud_frac_1,                               // In
                                  Real2D cloud_frac_2,                               // In
                                  Real2D mixt_frac,
                                  const precipitation_fractions& precip_fracs,         // In
                                  int cluster_allocation_strategy,                     // In
                                  bool l_lh_importance_sampling,                       // In
                                  bool l_lh_straight_mc,                               // In
                                  bool l_lh_clustered_sampling,                        // In
                                  bool l_lh_limit_weights,                             // In
                                  bool l_lh_var_frac,                                  // In
                                  bool l_lh_normalize_weights,                         // In
                                  bool l_lh_deterministic_test,                        // In
                                  bool l_calc_weights_all_levs_itime,                  // In
                                  Real4D X_u_all_levs, Real3D lh_sample_point_weights, // Out
                                  latin_hypercube_arrays_state& sampling_state,        // Case state
                                  mt95& rng,                                           // Case state
                                  const parameters_silhs_state& silhs                  // Case state
);
void compute_k_lh_start(ScratchWorkspace& workspace, const grid& gr, int nzt, int ngrdcol,
                        Real2D rcm_pdf,                                    // In
                        const pdf_parameter& pdf_params,                   // In
                        bool l_rcm_in_cloud_k_lh_start,                    // In
                        bool l_random_k_lh_start,                          // In
                        int lh_seed,                                       // In
                        Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // Out
                        mt95& rng                                          // Case state
);
void assert_consistent_cloud_frac(real_t chi_1, real_t chi_2,               // In
                                  real_t cloud_frac_1, real_t cloud_frac_2, // In
                                  real_t stdev_chi_1, real_t stdev_chi_2,   // In
                                  bool& l_error                             // Out
);
void assert_consistent_cf_component(real_t mu_chi_i, real_t sigma_chi_i, real_t cloud_frac_i, // In
                                    bool& l_error                                             // Out
);
void assert_correct_cloud_normal(int num_samples,            // In
                                 SilhsHostReal1D X_u_chi,    // In
                                 SilhsHostReal1D X_nl_chi,   // In
                                 SilhsHostInt1D X_mixt_comp, // In
                                 real_t cloud_frac_1,        // In
                                 real_t cloud_frac_2,        // In
                                 bool& l_error               // Out
);
void compute_arb_overlap(int nzt, int ngrdcol, int num_samples, int pdf_dim,
                         int d_uniform_extra,                                                 // In
                         Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, Real2D vert_corr, // In
                         Real4D rand_pool,                                                    // In
                         Real4D X_u_all_levs // InOut
);
void stats_accumulate_uniform_lh(
    ScratchWorkspace& workspace,                                                       // InOut
    int nzt, int num_samples, int ngrdcol,                                             // In
    Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft> l_in_precip_all_levs,            // In
    Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,                     // In
    Real3D X_u_chi_all_levs, const pdf_parameter& pdf_params,                          // In
    Real3D lh_sample_point_weights, Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // In
    stats_type& stats                                                                  // InOut
);
} // namespace

//-----------------------------------------------------------------------------
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
) {
  // intent(in)
  // intent(in)
  // intent(in)
  // intent(in)
  //rho_ds_zt, &
  // intent(in)
  // intent(in)
  // intent(in)
  // intent(in)
  // intent(inout)
  // intent(out)
  // intent(out)
  // intent(inout)
  // Description:
  //   Generate sample points of moisture, temperature, et cetera for the purpose
  //   of computing tendencies with a microphysics or radiation scheme.
  // References:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_silhs
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
  auto workspace_scope          = workspace.scope();
  constexpr int d_uniform_extra = 2;
  // dp1 selects PDF mixture component; dp2 selects precipitation membership.
  // The Native OpenACC-only straight-MC restriction is not needed by the
  // portable C++ host RNG plus Kokkos transform path; no CURAND state is used.
  if (silhs_config_flags.l_lh_importance_sampling && sequence_length != 1) {
    std::cerr << "Cloud weighted sampling requires sequence length be equal to 1.\n";
    throw std::runtime_error("Fatal error.");
  }
  if (silhs_config_flags.l_Lscale_vert_avg)
    throw std::runtime_error("CLUBB ERROR: l_Lscale_vert_avg has been depricated");
  bool l_error = false;
  Real4D X_u_all_levs =
      workspace_scope.view<Real4D>(Kokkos::WithoutInitializing, "X_u_all_levs", ngrdcol,
                                   num_samples, nzt, pdf_dim + d_uniform_extra);
  Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(Kokkos::WithoutInitializing,
                                                                   "k_lh_start", ngrdcol);
  Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft> l_in_precip =
      workspace_scope.view<Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "l_in_precip", ngrdcol, num_samples, nzt);
  Real4D Sigma_Cholesky1 = workspace_scope.view<Real4D>(
      Kokkos::WithoutInitializing, "Sigma_Cholesky1", pdf_dim, ngrdcol, nzt, pdf_dim);
  Real4D Sigma_Cholesky2 = workspace_scope.view<Real4D>(
      Kokkos::WithoutInitializing, "Sigma_Cholesky2", pdf_dim, ngrdcol, nzt, pdf_dim);
  Real3D cloud_frac = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing, "cloud_frac",
                                                   ngrdcol, num_samples, nzt);
  Real2D rcm_pdf =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "rcm_pdf", ngrdcol, nzt);
  Real2D X_vert_corr =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "X_vert_corr", ngrdcol, nzt);
  Real4D rand_pool = workspace_scope.view<Real4D>(Kokkos::WithoutInitializing, "rand_pool", ngrdcol,
                                                  num_samples, nzt, pdf_dim + d_uniform_extra);
  const auto rc_1  = pdf_params.rc_1;
  const auto rc_2  = pdf_params.rc_2;
  const auto mixt_frac     = pdf_params.mixt_frac;
  const auto cloud_frac_1  = pdf_params.cloud_frac_1;
  const auto cloud_frac_2  = pdf_params.cloud_frac_2;
  const auto precip_frac_1 = precip_fracs.precip_frac_1;
  const auto precip_frac_2 = precip_fracs.precip_frac_2;
  Kokkos::parallel_for(
      "generate_silhs_sample::rcm_pdf", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), k = (cell / ngrdcol);
        rcm_pdf(i, k) = compute_mean_binormal(rc_1(i, k), rc_2(i, k), mixt_frac(i, k));
      });
  compute_k_lh_start(workspace, gr, nzt, ngrdcol, rcm_pdf, pdf_params, // In
                     silhs_config_flags.l_rcm_in_cloud_k_lh_start,     // In
                     silhs_config_flags.l_random_k_lh_start,           // In
                     lh_seed,                                          // In
                     k_lh_start,                                       // Out
                     rng                                               // Case state
  );
  Kokkos::parallel_for(
      "generate_silhs_sample::Sigma", (ngrdcol * nzt * pdf_dim * pdf_dim), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), k = ((cell / ngrdcol) % nzt);
        const int j = ((cell / (ngrdcol * nzt)) % pdf_dim), p = (cell / (ngrdcol * nzt * pdf_dim));
        Sigma_Cholesky1(j, i, k, p) = corr_cholesky_mtx_1(i, k, p, j) * sigma1(i, k, p);
        Sigma_Cholesky2(j, i, k, p) = corr_cholesky_mtx_2(i, k, p, j) * sigma2(i, k, p);
      });
  const real_t grid_dir = gr.grid_dir;
  if (silhs_config_flags.l_max_overlap_in_cloud) {
    Kokkos::parallel_for(
        "generate_silhs_sample::cloud_overlap", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          if (rcm_pdf(i, k) > rc_tol)
            X_vert_corr(i, k) = one;
          else
            X_vert_corr(i, k) =
                Kokkos::exp(-vert_decorr_coef * (grid_dir * delta_zm(i, k) / Lscale(i, k)));
        });
  } else {
    Kokkos::parallel_for(
        "generate_silhs_sample::overlap", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          X_vert_corr(i, k) =
              Kokkos::exp(-vert_decorr_coef * (grid_dir * delta_zm(i, k) / Lscale(i, k)));
        });
  }
  if (clubb_at_least_debug_level_api(1)) {
    const auto host = to_host_space(workspace, X_vert_corr);
    bool invalid    = false;
    for (int k = 0; k < nzt; ++k)
      for (int i = 0; i < ngrdcol; ++i)
        invalid = invalid || host(i, k) > one || host(i, k) < zero;
    if (invalid) {
      std::cerr
          << "The vertical correlation in latin_hypercube_driveris not in the correct range\n";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          std::cerr << "i = " << (i + 1) << "k = " << (k + 1)
                    << "Vert. correlation = " << host(i, k) << '\n';
      throw std::runtime_error("Fatal error in vertical_overlap_driver");
    }
  }
  generate_random_pool(workspace, nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra, // In
                       lh_seed, gr,                                                    // In
                       silhs_config_flags.l_lh_deterministic_test,                     // In
                       rand_pool,                                                      // Out
                       rng                                                             // Case state
  );
  generate_all_uniform_samples(workspace, iter, pdf_dim, d_uniform_extra, num_samples,
                               sequence_length,                                  // In
                               nzt, ngrdcol, k_lh_start, X_vert_corr, rand_pool, // In
                               hm_metadata.iiPDF_chi,                            // In
                               cloud_frac_1,                                     // In
                               cloud_frac_2,                                     // In
                               mixt_frac, precip_fracs,                          // In
                               silhs_config_flags.cluster_allocation_strategy,   // In
                               silhs_config_flags.l_lh_importance_sampling,      // In
                               silhs_config_flags.l_lh_straight_mc,              // In
                               silhs_config_flags.l_lh_clustered_sampling,       // In
                               silhs_config_flags.l_lh_limit_weights,            // In
                               silhs_config_flags.l_lh_var_frac,                 // In
                               silhs_config_flags.l_lh_normalize_weights,        // In
                               silhs_config_flags.l_lh_deterministic_test,       // In
                               l_calc_weights_all_levs_itime,                    // In
                               X_u_all_levs, lh_sample_point_weights,            // Out
                               sampling_state, rng, silhs                        // Case state
  );
  Kokkos::parallel_for(
      "generate_silhs_sample::membership", (ngrdcol * num_samples * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                  k = (cell / (ngrdcol * num_samples));
        if (X_u_all_levs(i, sample, k, pdf_dim) < mixt_frac(i, k)) {
          X_mixt_comp_all_levs(i, sample, k) = 1;
          cloud_frac(i, sample, k)           = cloud_frac_1(i, k);
          if (X_u_all_levs(i, sample, k, (pdf_dim + 1)) < precip_frac_1(i, k))
            l_in_precip(i, sample, k) = true;
          else
            l_in_precip(i, sample, k) = false;
        } else {
          X_mixt_comp_all_levs(i, sample, k) = 2;
          cloud_frac(i, sample, k)           = cloud_frac_2(i, k);
          if (X_u_all_levs(i, sample, k, (pdf_dim + 1)) < precip_frac_2(i, k))
            l_in_precip(i, sample, k) = true;
          else
            l_in_precip(i, sample, k) = false;
        }
      });
  transform_uniform_samples_to_pdf(workspace, nzt, ngrdcol, num_samples, pdf_dim,
                                   d_uniform_extra,                  // In
                                   hm_metadata,                      // In
                                   Sigma_Cholesky1, Sigma_Cholesky2, // In
                                   mu1, mu2, X_mixt_comp_all_levs,   // In
                                   X_u_all_levs, cloud_frac,         // In
                                   l_in_precip,                      // In
                                   X_nl_all_levs                     // Out
  );
  if (stats && stats->l_sample) {
    stats_accumulate_uniform_lh(
        workspace, nzt, num_samples, ngrdcol, l_in_precip, // In
        X_mixt_comp_all_levs,
        to_execution_space(workspace, Kokkos::subview(X_u_all_levs, Kokkos::ALL, Kokkos::ALL,
                                                      Kokkos::ALL, hm_metadata.iiPDF_chi)),
        pdf_params,                          // In
        lh_sample_point_weights, k_lh_start, // In
        *stats                               // InOut
    );
  }
  if (stats && stats->l_last_sample) {
    if (l_output_2D_lognormal_dist)
      stats_lh_samples_write_lognormal(X_nl_all_levs, *stats, err_info);
    if (l_output_2D_uniform_dist)
      stats_lh_samples_write_uniform(X_u_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights,
                                     *stats, err_info);
    if (clubb_at_least_debug_level_api(0) && any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n'
                << "Fatal error writing SILHS 2D samples with stats in CLUBB SILHS procedure "
                   "generate_silhs_sample\n";
      return;
    }
  }
  if (clubb_at_least_debug_level_api(2)) {
    const auto X_u    = to_host_space(workspace, X_u_all_levs);
    const auto X_nl   = to_host_space(workspace, X_nl_all_levs);
    const auto X_mixt = to_host_space(workspace, X_mixt_comp_all_levs);
    bool invalid      = false;
    for (int p = 0; p < d_uniform_extra + pdf_dim; ++p)
      for (int k = 0; k < nzt; ++k)
        for (int sample = 0; sample < num_samples; ++sample)
          for (int i = 0; i < ngrdcol; ++i)
            invalid = invalid || X_u(i, sample, k, p) <= zero || X_u(i, sample, k, p) >= one;
    if (invalid) {
      std::cerr << "A uniform variate was not in the correct range.\n";
      l_error = true;
    }
    const auto chi_1   = to_host_space(workspace, pdf_params.chi_1),
               chi_2   = to_host_space(workspace, pdf_params.chi_2);
    const auto cf_1    = to_host_space(workspace, cloud_frac_1),
               cf_2    = to_host_space(workspace, cloud_frac_2);
    const auto stdev_1 = to_host_space(workspace, pdf_params.stdev_chi_1),
               stdev_2 = to_host_space(workspace, pdf_params.stdev_chi_2);
    for (int i = 0; i < ngrdcol; ++i) {
      for (int k = 0; k < nzt; ++k) {
        bool l_error_in_sub;
        assert_consistent_cloud_frac(chi_1(i, k), chi_2(i, k),     // In
                                     cf_1(i, k), cf_2(i, k),       // In
                                     stdev_1(i, k), stdev_2(i, k), // In
                                     l_error_in_sub                // Out
        );
        l_error = l_error || l_error_in_sub;
        assert_correct_cloud_normal(
            num_samples, // In
            to_host_space(workspace,
                          Kokkos::subview(X_u, i, Kokkos::ALL, k, hm_metadata.iiPDF_chi)), // In
            to_host_space(workspace,
                          Kokkos::subview(X_nl, i, Kokkos::ALL, k, hm_metadata.iiPDF_chi)), // In
            to_host_space(workspace, Kokkos::subview(X_mixt, i, Kokkos::ALL, k)),           // In
            cf_1(i, k),                                                                     // In
            cf_2(i, k),                                                                     // In
            l_error_in_sub                                                                  // Out
        );
        l_error = l_error || l_error_in_sub;
      }
    }
  }
  if (l_error)
    throw std::runtime_error("Fatal error in generate_silhs_sample");
}

namespace {

//-----------------------------------------------------------------------------
void generate_random_pool(ScratchWorkspace& workspace, int nzt, int ngrdcol, int pdf_dim,
                          int num_samples,
                          int d_uniform_extra,          // In
                          int lh_seed, const grid& gr,  // In
                          bool l_lh_deterministic_test, // In
                          Real4D rand_pool,             // Out
                          mt95& rng                     // Case state
) {
  // Description:
  //     This subroutine populates rand_pool with random numbers. There are
  //     different requirements for generating random numbers on a CPU vs a
  //     GPU, so the operation of this procedure depends on where or not
  //     we compile with -DCLUBB_GPU has been specified at compile time.
  //     We use a CUDA routine to generate on a GPU when available, but if
  //     -DCUDA is not specified then we generate randoms on the GPU and
  //     copy the results to the GPU.
  // References:
  //   clubb ticket 869
  // Author: Gunther Huebler
  // nzt: Number of vertical levels
  // ngrdcol: Number of grid columns
  // pdf_dim: Variates
  // num_samples: Number of samples
  // d_uniform_extra: Uniform variates included in uniform sample but not in normal/lognormal sample
  // gr: Grid variable type
  // l_lh_deterministic_test: Use a repeating overlap pool for repeatable testing
  // rand_pool: Pool of random reals, generated from integers
  const Kokkos::Array<real_t, 4> test_uniform_draws = {real_t(.125), real_t(.625), real_t(.375),
                                                       real_t(.875)};
  if (l_lh_deterministic_test) {
    Kokkos::parallel_for(
        "generate_random_pool::deterministic",
        (d_uniform_extra * ngrdcol * num_samples * nzt + ngrdcol * num_samples * nzt * pdf_dim),
        KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                    k                = ((cell / (ngrdcol * num_samples)) % nzt),
                    p                = (cell / (ngrdcol * num_samples * nzt));
          rand_pool(i, sample, k, p) = test_uniform_draws[((i + k + p + sample) % 4)];
        });
    return;
  }
  // Reinitialize every timestep to preserve Native restart reproducibility.
  // CPU MT95 remains the portable stream; Native CUDA CURAND is not delegated.
  rng.genrand_init_api(genrand_intg(lh_seed));
  auto workspace_scope = workspace.scope();
  auto host            = create_host_mirror(workspace, rand_pool);
  for (int p = 0; p < pdf_dim + d_uniform_extra; ++p)
    for (int k = gr.k_lb_zt; k * gr.grid_dir_indx <= gr.k_ub_zt * gr.grid_dir_indx;
         k += gr.grid_dir_indx)
      for (int sample = 0; sample < num_samples; ++sample)
        for (int i = 0; i < ngrdcol; ++i)
          host(i, sample, k, p) = rand_uniform_real(rng // Case state
          );
  Kokkos::deep_copy(rand_pool, host);
#ifdef SILHS_MULTI_COL_RAND_DUPLICATE
  Kokkos::parallel_for(
      "generate_random_pool::duplicate", ngrdcol * num_samples * nzt * (pdf_dim + d_uniform_extra),
      KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, sample = (cell / ngrdcol) % num_samples,
                  k = (cell / (ngrdcol * num_samples)) % nzt,
                  p = cell / (ngrdcol * num_samples * nzt);
        if (i >= 1)
          rand_pool(i, sample, k, p) = rand_pool(0, sample, k, p);
      });
#endif
}
} // namespace

namespace {

//-----------------------------------------------------------------------------
void generate_all_uniform_samples(ScratchWorkspace& workspace, int iter, int pdf_dim,
                                  int d_uniform_extra, int num_samples,
                                  int sequence_length, // In
                                  int nzt, int ngrdcol,
                                  Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // In
                                  Real2D X_vert_corr, Real4D rand_pool,              // In
                                  int iiPDF_chi,                                     // In
                                  Real2D cloud_frac_1,                               // In
                                  Real2D cloud_frac_2,                               // In
                                  Real2D mixt_frac,
                                  const precipitation_fractions& precip_fracs,         // In
                                  int cluster_allocation_strategy,                     // In
                                  bool l_lh_importance_sampling,                       // In
                                  bool l_lh_straight_mc,                               // In
                                  bool l_lh_clustered_sampling,                        // In
                                  bool l_lh_limit_weights,                             // In
                                  bool l_lh_var_frac,                                  // In
                                  bool l_lh_normalize_weights,                         // In
                                  bool l_lh_deterministic_test,                        // In
                                  bool l_calc_weights_all_levs_itime,                  // In
                                  Real4D X_u_all_levs, Real3D lh_sample_point_weights, // Out
                                  latin_hypercube_arrays_state& sampling_state,        // Case state
                                  mt95& rng,                                           // Case state
                                  const parameters_silhs_state& silhs                  // Case state
) {
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(in)
  // Intent(out)
  // Description:
  //   Generates uniform samples for all vertical levels, samples, and variates.
  //   Applys Latin Hypercude and importance sampling where conigured
  // References:
  //   V. E. Larson and D. P. Schanen, 2013. The Subgrid Importance Latin
  //   Hypercube Sampler (SILHS): a multivariate subcolumn generator
  // iter: Model iteration number
  // pdf_dim: Number of variates in CLUBB's PDF
  // d_uniform_extra: Uniform variates included in uniform sample but not in normal/lognormal sample
  // num_samples: Number of SILHS sample points
  // sequence_length: Number of timesteps before new sample points are picked
  // rand_pool: Array of randomly generated numbers
  // cloud_frac_1: The PDF parameters at k_lh_start
  // cloud_frac_2: The PDF parameters at k_lh_start
  // precip_fracs: Precipitation fractions      [-]
  // cluster_allocation_strategy: Strategy for distributing sample points
  // l_lh_importance_sampling: Do importance sampling (SILHS)
  // l_lh_straight_mc: Do not apply LH or importance sampling at all (SILHS)
  // l_lh_clustered_sampling: Use prescribed probability sampling with clusters (SILHS)
  // l_lh_limit_weights: Ensure weights stay under a given value
  // l_lh_var_frac: Prescribe variance fractions
  // l_lh_normalize_weights: Normalize weights to sum to num_samples
  // l_lh_deterministic_test: Ordered strata/repeating draws for repeatable testing
  // X_u_all_levs: Uniform sample at k_lh_start
  // lh_sample_point_weights: Weight of each sample point (all equal to one if importance sampling
  //   is not used)
  // Capture a local value for reference-taking Kokkos::max in the device paths.
  auto workspace_scope                   = workspace.scope();
  const real_t single_prec_thresh        = clubb_kokkos::single_prec_thresh;
  constexpr bool l_lh_old_cloud_weighted = false;
  if (l_lh_old_cloud_weighted && (num_samples % 2) != 0) {
    std::cerr << "Old cloud weighted sampling requires num_samples to be divisible by 2.\n";
    throw std::runtime_error("Fatal error.");
  }
  // Straight MC can populate all independent cells directly on the execution
  // space. Ordered LH/importance work shares case RNG/permutation on the host.
  if (!l_calc_weights_all_levs_itime) {
    if (l_lh_straight_mc) {
      Kokkos::parallel_for(
          "generate_all_uniform::mc",
          (d_uniform_extra * ngrdcol * num_samples * nzt + ngrdcol * num_samples * nzt * pdf_dim),
          KOKKOS_LAMBDA(int cell) {
            const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                      k = ((cell / (ngrdcol * num_samples)) % nzt),
                      p = (cell / (ngrdcol * num_samples * nzt));
            X_u_all_levs(i, sample, k, p) =
                Kokkos::max(single_prec_thresh,
                            Kokkos::min(one - single_prec_thresh, rand_pool(i, sample, k, p)));
          });
      Kokkos::deep_copy(lh_sample_point_weights, one);
    } else {
      const auto start = to_host_space(workspace, k_lh_start);
      const auto cf1   = to_host_space(workspace, cloud_frac_1),
                 cf2   = to_host_space(workspace, cloud_frac_2),
                 mf    = to_host_space(workspace, mixt_frac);
      const auto pf1   = to_host_space(workspace, precip_fracs.precip_frac_1),
                 pf2   = to_host_space(workspace, precip_fracs.precip_frac_2);
      auto X           = create_host_mirror(workspace, X_u_all_levs);
      auto weights     = create_host_mirror(workspace, lh_sample_point_weights);
      for (int i = 0; i < ngrdcol; ++i) {
        const int k           = start(i);
        SilhsHostReal2D X_one = workspace_scope.view<SilhsHostReal2D>(
            Kokkos::WithoutInitializing, "X_u_one_lev", num_samples, pdf_dim + d_uniform_extra);
        generate_uniform_lh_sample(iter, num_samples, sequence_length,
                                   (d_uniform_extra + pdf_dim), // In
                                   l_lh_deterministic_test,     // In
                                   X_one,                       // Out
                                   sampling_state, rng          // Case state
        );
        SilhsHostReal1D weight_one = workspace_scope.view<SilhsHostReal1D>(
            Kokkos::WithoutInitializing, "weights_one_lev", num_samples);
        if (l_lh_importance_sampling) {
          if (l_lh_old_cloud_weighted) {
            const auto matrix = to_host_space(workspace, sampling_state.one_height_time_matrix);
            auto chi = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi));
            auto dp1 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim));
            cloud_weighted_sampling_driver(
                workspace, num_samples,
                to_host_space(workspace, Kokkos::subview(matrix, Kokkos::ALL, iiPDF_chi)),
                to_host_space(workspace, Kokkos::subview(matrix, Kokkos::ALL, pdf_dim)), // In
                cf1(i, k), cf2(i, k), mf(i, k),                                          // In
                chi, dp1,                                                                // InOut
                weight_one,                                                              // Out
                rng // Case state
            );
            Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi), chi);
            Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, pdf_dim), dp1);
          } else {
            auto chi = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi));
            auto dp1 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim));
            auto dp2 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim + 1));
            importance_sampling_driver(workspace, num_samples,                               // In
                                       cf1(i, k), cf2(i, k),                                 // In
                                       mf(i, k),                                             // In
                                       pf1(i, k), pf2(i, k),                                 // In
                                       cluster_allocation_strategy, l_lh_clustered_sampling, // In
                                       l_lh_limit_weights, l_lh_var_frac,
                                       l_lh_normalize_weights, // In
                                       chi, dp1, dp2,          // InOut
                                       weight_one,             // Out
                                       rng, silhs              // Case state
            );
            Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi), chi);
            Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, pdf_dim), dp1);
            Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, (pdf_dim + 1)), dp2);
          }
          // Native clips only this importance branch before arbitrary overlap.
          for (int p = 0; p < d_uniform_extra + pdf_dim; ++p)
            for (int sample = 0; sample < num_samples; ++sample)
              X_one(sample, p) = std::max(single_prec_thresh,
                                          std::min(one - single_prec_thresh, X_one(sample, p)));
          for (int level = 0; level < nzt; ++level)
            for (int sample = 0; sample < num_samples; ++sample)
              weights(i, sample, level) = weight_one(sample);
        } else {
          for (int level = 0; level < nzt; ++level)
            for (int sample = 0; sample < num_samples; ++sample)
              weights(i, sample, level) = one;
        }
        for (int p = 0; p < d_uniform_extra + pdf_dim; ++p)
          for (int sample = 0; sample < num_samples; ++sample)
            X(i, sample, k, p) = X_one(sample, p);
      }
      Kokkos::deep_copy(X_u_all_levs, X);
      Kokkos::deep_copy(lh_sample_point_weights, weights);
    }
    compute_arb_overlap(nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra, // In
                        k_lh_start, X_vert_corr, rand_pool,                  // In
                        X_u_all_levs                                         // InOut
    );
  } else {
    if (l_lh_straight_mc) {
      Kokkos::parallel_for(
          "generate_all_uniform::independent_mc",
          (d_uniform_extra * ngrdcol * num_samples * nzt + ngrdcol * num_samples * nzt * pdf_dim),
          KOKKOS_LAMBDA(int cell) {
            const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                      k = ((cell / (ngrdcol * num_samples)) % nzt),
                      p = (cell / (ngrdcol * num_samples * nzt));
            X_u_all_levs(i, sample, k, p) =
                Kokkos::max(single_prec_thresh,
                            Kokkos::min(one - single_prec_thresh, rand_pool(i, sample, k, p)));
          });
      Kokkos::deep_copy(lh_sample_point_weights, one);
    } else {
      const auto cf1 = to_host_space(workspace, cloud_frac_1),
                 cf2 = to_host_space(workspace, cloud_frac_2),
                 mf  = to_host_space(workspace, mixt_frac);
      const auto pf1 = to_host_space(workspace, precip_fracs.precip_frac_1),
                 pf2 = to_host_space(workspace, precip_fracs.precip_frac_2);
      auto X         = create_host_mirror(workspace, X_u_all_levs);
      auto weights   = create_host_mirror(workspace, lh_sample_point_weights);
      for (int i = 0; i < ngrdcol; ++i) {
        for (int k = 0; k < nzt; ++k) {
          SilhsHostReal2D X_one = workspace_scope.view<SilhsHostReal2D>(
              Kokkos::WithoutInitializing, "X_u_one_lev", num_samples, pdf_dim + d_uniform_extra);
          generate_uniform_lh_sample(iter, num_samples, sequence_length,
                                     (d_uniform_extra + pdf_dim), // In
                                     l_lh_deterministic_test,     // In
                                     X_one,                       // Out
                                     sampling_state, rng          // Case state
          );
          SilhsHostReal1D weight_one = workspace_scope.view<SilhsHostReal1D>(
              Kokkos::WithoutInitializing, "weights_one_lev", num_samples);
          if (l_lh_importance_sampling) {
            if (l_lh_old_cloud_weighted) {
              const auto matrix = to_host_space(workspace, sampling_state.one_height_time_matrix);
              auto chi = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi));
              auto dp1 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim));
              cloud_weighted_sampling_driver(
                  workspace, num_samples,
                  to_host_space(workspace, Kokkos::subview(matrix, Kokkos::ALL, iiPDF_chi)),
                  to_host_space(workspace, Kokkos::subview(matrix, Kokkos::ALL, pdf_dim)), // In
                  cf1(i, k), cf2(i, k), mf(i, k),                                          // In
                  chi, dp1,                                                                // InOut
                  weight_one,                                                              // Out
                  rng // Case state
              );
              Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi), chi);
              Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, pdf_dim), dp1);
            } else {
              auto chi = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi));
              auto dp1 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim));
              auto dp2 = to_host_space(workspace, Kokkos::subview(X_one, Kokkos::ALL, pdf_dim + 1));
              importance_sampling_driver(workspace, num_samples,                               // In
                                         cf1(i, k), cf2(i, k),                                 // In
                                         mf(i, k),                                             // In
                                         pf1(i, k), pf2(i, k),                                 // In
                                         cluster_allocation_strategy, l_lh_clustered_sampling, // In
                                         l_lh_limit_weights, l_lh_var_frac,
                                         l_lh_normalize_weights, // In
                                         chi, dp1, dp2,          // InOut
                                         weight_one,             // Out
                                         rng, silhs              // Case state
              );
              Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, iiPDF_chi), chi);
              Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, pdf_dim), dp1);
              Kokkos::deep_copy(Kokkos::subview(X_one, Kokkos::ALL, (pdf_dim + 1)), dp2);
            }
          } else
            Kokkos::deep_copy(weight_one, one);
          // This branch clips even when importance sampling is disabled.
          for (int p = 0; p < d_uniform_extra + pdf_dim; ++p)
            for (int sample = 0; sample < num_samples; ++sample)
              X(i, sample, k, p) = std::max(single_prec_thresh,
                                            std::min(one - single_prec_thresh, X_one(sample, p)));
          for (int sample = 0; sample < num_samples; ++sample)
            weights(i, sample, k) = weight_one(sample);
        }
      }
      Kokkos::deep_copy(X_u_all_levs, X);
      Kokkos::deep_copy(lh_sample_point_weights, weights);
    }
  }
}
} // namespace

namespace {

//-----------------------------------------------------------------------------
void compute_k_lh_start(ScratchWorkspace& workspace, const grid& gr, int nzt, int ngrdcol,
                        Real2D rcm_pdf,                                    // In
                        const pdf_parameter& pdf_params,                   // In
                        bool l_rcm_in_cloud_k_lh_start,                    // In
                        bool l_random_k_lh_start,                          // In
                        int lh_seed,                                       // In
                        Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // Out
                        mt95& rng                                          // Case state
) {
  // Description:
  //   Determines the starting SILHS sample level
  // References:
  //   None
  // gr: Grid variable type
  // nzt: Number of vertical levels
  // ngrdcol: Number of grid columns
  // rcm_pdf: Liquid water mixing ratio               [kg/kg]
  // pdf_params: PDF parameters       [units vary]
  // l_rcm_in_cloud_k_lh_start: Determine k_lh_start based on maximum within-cloud rcm
  // l_random_k_lh_start: k_lh_start found randomly between max rcm and rcm_in_cloud
  // lh_seed: Random number generator seed
  // k_lh_start: Zero-based starting SILHS sample slot
  // Kokkos::max needs an addressable device copy of this Native threshold.
  auto workspace_scope        = workspace.scope();
  const real_t cloud_frac_min = clubb_kokkos::cloud_frac_min;
  Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start_rcm_in_cloud =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "k_lh_start_rcm_in_cloud", ngrdcol);
  Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start_rcm =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(Kokkos::WithoutInitializing,
                                                                   "k_lh_start_rcm", ngrdcol);
  Real2D cloud_frac_pdf =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "cloud_frac_pdf", ngrdcol, nzt);
  const int grid_dir_indx = gr.grid_dir_indx;
  if (l_rcm_in_cloud_k_lh_start || l_random_k_lh_start) {
    const auto cf1 = pdf_params.cloud_frac_1, cf2 = pdf_params.cloud_frac_2,
               mf = pdf_params.mixt_frac;
    Kokkos::parallel_for(
        "compute_k_lh_start::cloud_frac", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          cloud_frac_pdf(i, k) = compute_mean_binormal(cf1(i, k), cf2(i, k), mf(i, k));
        });
    Kokkos::parallel_for(
        "compute_k_lh_start::in_cloud", ngrdcol, KOKKOS_LAMBDA(int i) {
          k_lh_start_rcm_in_cloud(i) = nzt / 2 - 1;
          if (grid_dir_indx < 0)
            k_lh_start_rcm_in_cloud(i) = nzt - k_lh_start_rcm_in_cloud(i) - 1;
          real_t rcm_in_cloud_max = zero;
          for (int k = 0; k < nzt; ++k) {
            const real_t rcm_in_cloud =
                rcm_pdf(i, k) / Kokkos::max(cloud_frac_pdf(i, k), cloud_frac_min);
            if (rcm_in_cloud > rcm_in_cloud_max) {
              rcm_in_cloud_max           = rcm_in_cloud;
              k_lh_start_rcm_in_cloud(i) = k;
            }
          }
        });
  }
  if (!l_rcm_in_cloud_k_lh_start || l_random_k_lh_start) {
    Kokkos::parallel_for(
        "compute_k_lh_start::gridbox", ngrdcol, KOKKOS_LAMBDA(int i) {
          k_lh_start_rcm(i) = nzt / 2 - 1;
          if (grid_dir_indx < 0)
            k_lh_start_rcm(i) = nzt - k_lh_start_rcm(i) - 1;
          real_t rcm_max = zero;
          for (int k = 0; k < nzt; ++k)
            if (rcm_pdf(i, k) > rcm_max) {
              rcm_max           = rcm_pdf(i, k);
              k_lh_start_rcm(i) = k;
            }
        });
  }
  if (l_random_k_lh_start) {
    const auto in_cloud = to_host_space(workspace, k_lh_start_rcm_in_cloud),
               gridbox  = to_host_space(workspace, k_lh_start_rcm);
    auto start          = create_host_mirror(workspace, k_lh_start);
    rng.genrand_init_api(genrand_intg(lh_seed));
    // Shifting both bounds preserves high-low+1 and the exact signed-MT95 draw.
    for (int i = 0; i < ngrdcol; ++i) {
      if (in_cloud(i) > gridbox(i))
        start(i) = rand_integer_in_range(gridbox(i), in_cloud(i), // In
                                         rng                      // Case state
        );
      else if (gridbox(i) > in_cloud(i))
        start(i) = rand_integer_in_range(in_cloud(i), gridbox(i), // In
                                         rng                      // Case state
        );
      else
        start(i) = gridbox(i); // Equal bounds consume no draw in Native.
      if (grid_dir_indx < 0)
        start(i) = in_cloud(i) + gridbox(i) - start(i);
    }
    Kokkos::deep_copy(k_lh_start, start);
  } else if (l_rcm_in_cloud_k_lh_start)
    Kokkos::deep_copy(k_lh_start, k_lh_start_rcm_in_cloud);
  else
    Kokkos::deep_copy(k_lh_start, k_lh_start_rcm);
}
} // namespace

//-----------------------------------------------------------------------------
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
) {
  // In
  // In
  // In
  // InOut
  // In
  // Out
  // Out
  // Out
  // Description:
  //   Derives from the SILHS sampling structure X_nl_all_levs the variables
  //   rt, thl, rc, rv, and Nc, for all sample points and height levels.
  // References:
  //   ticket:751
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
  // Capture local values for reference-taking Kokkos::max in the device kernel.
  auto workspace_scope                   = workspace.scope();
  const real_t rt_tol                    = clubb_kokkos::rt_tol;
  const real_t zero                      = clubb_kokkos::zero;
  constexpr bool l_clip_hydromet_samples = false;
  chi_eta_2_rtthl(
      nzt, ngrdcol, num_samples,           // In
      pdf_params.rt_1, pdf_params.thl_1,   // In
      pdf_params.rt_2, pdf_params.thl_2,   // In
      pdf_params.crt_1, pdf_params.cthl_1, // In
      pdf_params.crt_2, pdf_params.cthl_2, // In
      pdf_params.chi_1, pdf_params.chi_2,  // In
      to_execution_space(workspace, Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL,
                                                    Kokkos::ALL, hm_metadata.iiPDF_chi)),
      to_execution_space(workspace, Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL,
                                                    Kokkos::ALL, hm_metadata.iiPDF_eta)), // In
      X_mixt_comp_all_levs,                                                               // In
      lh_rt_clipped, lh_thl_clipped                                                       // Out
  );
  const int iiPDF_chi = hm_metadata.iiPDF_chi, iiPDF_Ncn = hm_metadata.iiPDF_Ncn;
  Kokkos::parallel_for(
      "clip_transform_silhs_output", (ngrdcol * num_samples * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                  k                 = (cell / (ngrdcol * num_samples));
        lh_rt_clipped(i, sample, k) = Kokkos::max(lh_rt_clipped(i, sample, k), rt_tol);
        lh_rc_clipped(i, sample, k) = Kokkos::max(X_nl_all_levs(i, sample, k, iiPDF_chi), zero);
        lh_rc_clipped(i, sample, k) =
            Kokkos::min(lh_rc_clipped(i, sample, k), lh_rt_clipped(i, sample, k) - rt_tol);
        lh_rv_clipped(i, sample, k) = lh_rt_clipped(i, sample, k) - lh_rc_clipped(i, sample, k);
        if (l_use_Ncn_to_Nc) {
          if (X_nl_all_levs(i, sample, k, iiPDF_chi) > zero)
            lh_Nc_clipped(i, sample, k) = X_nl_all_levs(i, sample, k, iiPDF_Ncn);
          else
            lh_Nc_clipped(i, sample, k) = zero;
        } else
          lh_Nc_clipped(i, sample, k) = X_nl_all_levs(i, sample, k, iiPDF_Ncn);
      });
  if (l_clip_hydromet_samples) {
    // Native local option is false. Retain its complete concentration-clipping
    // sequence through the existing permanent hydrometeor radius-limit API.
    Real2D hydromet_pts = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "hydromet_pts",
                                                       nzt, hydromet_dim);
    Real2D hydromet_pts_clipped = workspace_scope.view<Real2D>(
        Kokkos::WithoutInitializing, "hydromet_pts_clipped", nzt, hydromet_dim);
    for (int sample = 0; sample < num_samples; ++sample) {
      for (int i = 0; i < ngrdcol; ++i) {
        for (int hm_idx = 0; hm_idx < hydromet_dim; ++hm_idx) {
          const int pdf_idx = ((hydromet2pdf_idx(hm_idx, hm_metadata) + 1)) - 1;
          Kokkos::parallel_for(
              "clip_transform::pack", nzt, KOKKOS_LAMBDA(int k) {
                hydromet_pts(k, hm_idx) = X_nl_all_levs(i, sample, k, pdf_idx);
              });
        }
        clip_hydromet_conc_mvr(workspace, nzt, hydromet_dim, hm_metadata, hydromet_pts,
                               hydromet_pts_clipped);
        for (int hm_idx = 0; hm_idx < hydromet_dim; ++hm_idx) {
          const int pdf_idx = ((hydromet2pdf_idx(hm_idx, hm_metadata) + 1)) - 1;
          Kokkos::parallel_for(
              "clip_transform::unpack", nzt, KOKKOS_LAMBDA(int k) {
                X_nl_all_levs(i, sample, k, pdf_idx) = hydromet_pts_clipped(k, hm_idx);
              });
        }
      }
    }
  }
}

namespace {

//-----------------------------------------------------------------------------
void assert_consistent_cloud_frac(real_t chi_1, real_t chi_2,               // In
                                  real_t cloud_frac_1, real_t cloud_frac_2, // In
                                  real_t stdev_chi_1, real_t stdev_chi_2,   // In
                                  bool& l_error                             // Out
) {
  // Description:
  //   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  //   stdev_chi_i in pdf_params for each PDF component.
  // References:
  //   Eric Raut
  // chi_1: Discription of the distribution of chi
  // chi_2: Discription of the distribution of chi
  // cloud_frac_1: and cloud fraction [units vary]
  // cloud_frac_2: and cloud fraction [units vary]
  // l_error: True if the assertion check fails
  l_error = false;
  bool l_error_in_sub;
  assert_consistent_cf_component(chi_1, stdev_chi_1, cloud_frac_1, // In
                                 l_error_in_sub                    // Out
  );
  l_error = l_error || l_error_in_sub;
  if (l_error_in_sub)
    std::cerr << "Cloud fraction is inconsistent in PDF component 1\n";
  assert_consistent_cf_component(chi_2, stdev_chi_2, cloud_frac_2, // In
                                 l_error_in_sub                    // Out
  );
  l_error = l_error || l_error_in_sub;
  if (l_error_in_sub)
    std::cerr << "Cloud fraction is inconsistent in PDF component 2\n";
}
} // namespace

namespace {

//-----------------------------------------------------------------------------
void assert_consistent_cf_component(real_t mu_chi_i, real_t sigma_chi_i, real_t cloud_frac_i, // In
                                    bool& l_error                                             // Out
) {
  // Description:
  //   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  //   stdev_chi_i for a PDF component.
  //   The SILHS sample generation process relies on precisely a cloud_frac
  //   amount of mass in the cloudy portion of the PDF of chi, that is, where
  //   chi > 0. In other words, the probability that chi > 0 should be exactly
  //   cloud_frac.
  //   Stated even more mathematically, CDF_chi(0) = 1 - cloud_frac, where
  //   CDF_chi is the cumulative distribution function of chi. This can be
  //   expressed as invCDF_chi(1 - cloud_frac) = zero.
  //   This subroutine uses ltqnorm, which is apparently a fancy name for the
  //   inverse cumulative distribution function of the standard normal
  //   distribution.
  // References:
  //   Eric Raut
  // mu_chi_i: Mean of chi in a PDF component
  // sigma_chi_i: Standard deviation of chi in a PDF component
  // cloud_frac_i: Cloud fraction in a PDF component
  // l_error: True if the assertion check fails
  constexpr real_t ltqnorm_min_arg = real_t(1.e-5), box_half_width = real_t(5.e-6);
  l_error = false;
  // Keep this omission identical to closure's zero-cloud special condition.
  if (std::abs(mu_chi_i) <= eps && sigma_chi_i <= chi_tol)
    return;
  real_t one_minus_ltqnorm_arg = cloud_frac_i - box_half_width;
  if (one_minus_ltqnorm_arg >= ltqnorm_min_arg) {
    if (one_minus_ltqnorm_arg > one - ltqnorm_min_arg)
      one_minus_ltqnorm_arg = one - ltqnorm_min_arg;
    const real_t chi_std_normal = ltqnorm(one - one_minus_ltqnorm_arg // In
    );
    const real_t chi            = chi_std_normal * sigma_chi_i + mu_chi_i;
    if (chi <= zero) {
      l_error = true;
      std::cerr << "chi (left side of box) = " << chi << '\n';
    }
  }
  one_minus_ltqnorm_arg = cloud_frac_i + box_half_width;
  if (one_minus_ltqnorm_arg <= one - ltqnorm_min_arg) {
    if (one_minus_ltqnorm_arg < ltqnorm_min_arg)
      one_minus_ltqnorm_arg = ltqnorm_min_arg;
    const real_t chi_std_normal = ltqnorm(one - one_minus_ltqnorm_arg // In
    );
    const real_t chi            = chi_std_normal * sigma_chi_i + mu_chi_i;
    if (chi > zero) {
      l_error = true;
      std::cerr << "chi (right side of box) = " << chi << '\n';
    }
  }
  if (l_error)
    std::cerr << "In assert_consistent_cf_component, cloud_frac_i is inconsistent with mu_chi_i "
                 "and stdev_chi_i.\n"
              << "mu_chi_i = " << mu_chi_i << "\nsigma_chi_i = " << sigma_chi_i
              << "\ncloud_frac_i = " << cloud_frac_i << '\n';
}
} // namespace

namespace {

//-----------------------------------------------------------------------------
void assert_correct_cloud_normal(int num_samples,            // In
                                 SilhsHostReal1D X_u_chi,    // In
                                 SilhsHostReal1D X_nl_chi,   // In
                                 SilhsHostInt1D X_mixt_comp, // In
                                 real_t cloud_frac_1,        // In
                                 real_t cloud_frac_2,        // In
                                 bool& l_error               // Out
) {
  // Description:
  //   Asserts that all SILHS sample points that are in cloud in uniform space
  //   are in cloud in normal space, and that all SILHS sample points that are
  //   in clear air in uniform space are in clear air in normal space.
  // References:
  //   None
  // num_samples: Number of SILHS sample points
  // X_u_chi: Samples of chi in uniform space
  // X_nl_chi: Samples of chi in normal space
  // X_mixt_comp: PDF component of each sample
  // cloud_frac_1: Cloud fraction in PDF component 1
  // cloud_frac_2: Cloud fraction in PDF component 2
  // l_error: True if the assertion check fails
  l_error = false;
  // Valid caller component labels are 1 or2; Native leaves cloud_frac_i
  // unwritten for another label. No extra fallback or validation is invented.
  real_t cloud_frac_i;
  for (int sample = 0; sample < num_samples; ++sample) {
    if (X_mixt_comp(sample) == 1)
      cloud_frac_i = cloud_frac_1;
    else if (X_mixt_comp(sample) == 2)
      cloud_frac_i = cloud_frac_2;
    if (X_u_chi(sample) < one - cloud_frac_i) {
      if (X_nl_chi(sample) > real_t(1000) * single_prec_thresh) {
        l_error = true;
        std::cerr << "X_nl_chi(" << (sample + 1) << ") > " << single_prec_thresh << '\n';
      }
    } else if (X_u_chi(sample) >= one - cloud_frac_i && X_u_chi(sample) < one) {
      if (X_nl_chi(sample) <= -real_t(1000) * single_prec_thresh) {
        l_error = true;
        std::cerr << "X_nl_chi(" << (sample + 1) << ") <= " << -single_prec_thresh << '\n';
      }
    } else
      throw std::runtime_error("X_u_chi not in correct range in assert_correct_cloud_normal");
  }
  if (l_error) {
    std::cerr << "In assert_correct_cloud_normal:\n"
              << "The 'cloudiness' of points in uniform and normal space is not consistent\n"
              << "    X_u_chi         X_nl_chi \n";
    for (int sample = 0; sample < num_samples; ++sample)
      std::cerr << std::setw(4) << (sample + 1) << std::setw(20) << std::setprecision(4)
                << X_u_chi(sample) << std::setw(20) << X_nl_chi(sample) << '\n';
    l_error = true;
  }
}
} // namespace

//-----------------------------------------------------------------------------
void latin_hypercube_2D_output_api(int nzt,
                                   Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> zt,
                                   int pdf_dim, int num_samples,        // In
                                   const hm_metadata_type& hm_metadata, // In
                                   stats_type* stats,                   // Optional InOut
                                   err_info_type& err_info              // InOut
) {
  // err_info: err_info struct containing err_code and err_header
  const int iiPDF_chi = hm_metadata.iiPDF_chi, iiPDF_eta = hm_metadata.iiPDF_eta,
            iiPDF_w  = hm_metadata.iiPDF_w;
  const int iiPDF_rr = hm_metadata.iiPDF_rr, iiPDF_ri = hm_metadata.iiPDF_ri,
            iiPDF_rs = hm_metadata.iiPDF_rs, iiPDF_rg = hm_metadata.iiPDF_rg;
  const int iiPDF_Nr = hm_metadata.iiPDF_Nr, iiPDF_Ni = hm_metadata.iiPDF_Ni,
            iiPDF_Ns = hm_metadata.iiPDF_Ns, iiPDF_Ng = hm_metadata.iiPDF_Ng,
            iiPDF_Ncn = hm_metadata.iiPDF_Ncn;
  if (!stats)
    return;
  if (!stats->enabled)
    return;
  std::vector<std::string> nl_variable_names, u_variable_names;
  if (l_output_2D_lognormal_dist) {
    nl_variable_names.resize(pdf_dim);
    nl_variable_names[iiPDF_chi] = "chi";
    nl_variable_names[iiPDF_eta] = "eta";
    nl_variable_names[iiPDF_w]   = "w";
    if (iiPDF_rr > -1)
      nl_variable_names[iiPDF_rr] = "rr";
    if (iiPDF_ri > -1)
      nl_variable_names[iiPDF_ri] = "ri";
    if (iiPDF_rs > -1)
      nl_variable_names[iiPDF_rs] = "rs";
    if (iiPDF_rg > -1)
      nl_variable_names[iiPDF_rg] = "rg";
    if (iiPDF_Nr > -1)
      nl_variable_names[iiPDF_Nr] = "Nr";
    if (iiPDF_Ncn > -1)
      nl_variable_names[iiPDF_Ncn] = "Ncn";
    if (iiPDF_Ni > -1)
      nl_variable_names[iiPDF_Ni] = "Ni";
    if (iiPDF_Ns > -1)
      nl_variable_names[iiPDF_Ns] = "Ns";
    if (iiPDF_Ng > -1)
      nl_variable_names[iiPDF_Ng] = "Ng";
  } else
    nl_variable_names.resize(0);
  if (l_output_2D_uniform_dist) {
    u_variable_names.resize((pdf_dim + 4));
    u_variable_names[iiPDF_chi] = "chi";
    u_variable_names[iiPDF_eta] = "eta";
    u_variable_names[iiPDF_w]   = "w";
    if (iiPDF_rr > -1)
      u_variable_names[iiPDF_rr] = "rr";
    if (iiPDF_ri > -1)
      u_variable_names[iiPDF_ri] = "ri";
    if (iiPDF_rs > -1)
      u_variable_names[iiPDF_rs] = "rs";
    if (iiPDF_rg > -1)
      u_variable_names[iiPDF_rg] = "rg";
    if (iiPDF_Nr > -1)
      u_variable_names[iiPDF_Nr] = "Nr";
    if (iiPDF_Ncn > -1)
      u_variable_names[iiPDF_Ncn] = "Ncn";
    if (iiPDF_Ni > -1)
      u_variable_names[iiPDF_Ni] = "Ni";
    if (iiPDF_Ns > -1)
      u_variable_names[iiPDF_Ns] = "Ns";
    if (iiPDF_Ng > -1)
      u_variable_names[iiPDF_Ng] = "Ng";
    int p               = pdf_dim;
    u_variable_names[p] = "dp1";
    p                   = (pdf_dim + 1);
    u_variable_names[p] = "dp2";
    p                   = (pdf_dim + 2);
    u_variable_names[p] = "X_mixt_comp";
    p                   = (pdf_dim + 3);
    u_variable_names[p] = "lh_sample_point_weights";
  } else
    u_variable_names.resize(0);
  stats_lh_samples_init(num_samples, nzt, nl_variable_names, u_variable_names, zt, *stats,
                        err_info);
  // Stats finalization closes output. Only local name arrays are released here.
  nl_variable_names.clear();
  u_variable_names.clear();
  if (any_equal(err_info.err_code, clubb_fatal_error)) {
    std::cerr << err_info.err_header_global
              << "\nFatal error configuring SILHS 2D output in latin_hypercube_2D_output_api\n";
    return;
  }
}

namespace {

//-----------------------------------------------------------------------------
void compute_arb_overlap(int nzt, int ngrdcol, int num_samples, int pdf_dim,
                         int d_uniform_extra,                                                 // In
                         Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, Real2D vert_corr, // In
                         Real4D rand_pool,                                                    // In
                         Real4D X_u_all_levs // InOut
) {
  // Description:
  //   Re-computes X_u (uniform sample) using an arbitrary correlation specified
  //   by X_vert_corr (which can vary with height).
  //   This is an improved algorithm that doesn't require us to convert from a
  //   unifrom distribution to a Gaussian distribution and back again.
  // References:
  //   None

  // Each (column,sample,variate) retains two ordered vertical recurrences.
  // nzt: Vertical levels
  // ngrdcol: Columns
  // num_samples: Number of SILHS sample points
  // pdf_dim: Number of variates in CLUBB's PDF
  // d_uniform_extra: Uniform variates included in uniform sample but not in normal/lognormal sample
  // vert_corr: Vertical correlation between k points in range [0,1]   [-]
  // rand_pool: Array of randomly generated numbers
  // X_u_all_levs: Uniform distribution of 1 variate at all levels [-] The value of this variate at
  //   k_lh_start should already be populated in this array and will be used to fill in the other
  //   levels.
  Kokkos::parallel_for(
      "compute_arb_overlap::upward",
      (d_uniform_extra * ngrdcol * num_samples + ngrdcol * num_samples * pdf_dim),
      KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, sample = (cell / ngrdcol) % num_samples,
                  p            = cell / (ngrdcol * num_samples);
        const int k_lh_start_i = k_lh_start(i);
        real_t unbounded_point = X_u_all_levs(i, sample, k_lh_start_i, p);
        for (int k = k_lh_start_i; k < (nzt - 1); ++k) {
          const real_t half_width = one - vert_corr(i, (k + 1));
          const real_t min_val    = unbounded_point - half_width;
          const real_t offset     = two * half_width * rand_pool(i, sample, (k + 1), p);
          unbounded_point         = min_val + offset;
          if (unbounded_point > one - single_prec_thresh)
            unbounded_point = two - unbounded_point - two * single_prec_thresh;
          else if (unbounded_point < single_prec_thresh)
            unbounded_point = -unbounded_point + two * single_prec_thresh;
          X_u_all_levs(i, sample, (k + 1), p) = unbounded_point;
        }
      });
  Kokkos::parallel_for(
      "compute_arb_overlap::downward",
      (d_uniform_extra * ngrdcol * num_samples + ngrdcol * num_samples * pdf_dim),
      KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, sample = (cell / ngrdcol) % num_samples,
                  p            = cell / (ngrdcol * num_samples);
        const int k_lh_start_i = k_lh_start(i);
        real_t unbounded_point = X_u_all_levs(i, sample, k_lh_start_i, p);
        for (int k = k_lh_start_i; k >= 1; --k) {
          const real_t half_width = one - vert_corr(i, (k - 1));
          const real_t min_val    = unbounded_point - half_width;
          const real_t offset     = two * half_width * rand_pool(i, sample, (k - 1), p);
          unbounded_point         = min_val + offset;
          if (unbounded_point > one - single_prec_thresh)
            unbounded_point = two - unbounded_point - two * single_prec_thresh;
          else if (unbounded_point < single_prec_thresh)
            unbounded_point = -unbounded_point + two * single_prec_thresh;
          X_u_all_levs(i, sample, (k - 1), p) = unbounded_point;
        }
      });
}
} // namespace

//-----------------------------------------------------------------------------
void stats_accumulate_lh_api(ScratchWorkspace& workspace, const grid& gr, int nzt, int ngrdcol,
                             int num_samples, int pdf_dim,
                             Real2D rho_ds_zt,                                      // In
                             int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                             Real3D lh_sample_point_weights, Real4D X_nl_all_levs,  // In
                             Real3D lh_rt_clipped, Real3D lh_thl_clipped,           // In
                             Real3D lh_rc_clipped, Real3D lh_rv_clipped,            // In
                             Real3D lh_Nc_clipped,                                  // In
                             stats_type* stats // Optional InOut
) {
  // Description:
  //   Clip subcolumns from latin hypercube and create stats for diagnostic
  //   purposes.
  // References:
  //   None
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
  const int iirr = hm_metadata.iirr, iirs = hm_metadata.iirs, iiri = hm_metadata.iiri,
            iirg = hm_metadata.iirg;
  const int iiNr = hm_metadata.iiNr, iiNs = hm_metadata.iiNs, iiNi = hm_metadata.iiNi,
            iiNg      = hm_metadata.iiNg;
  const int iiPDF_chi = hm_metadata.iiPDF_chi, iiPDF_eta = hm_metadata.iiPDF_eta,
            iiPDF_w = hm_metadata.iiPDF_w, iiPDF_Ncn = hm_metadata.iiPDF_Ncn;
  bool l_stats = false;
  if (stats)
    l_stats = stats->l_sample;
  if (!l_stats)
    return;
  auto lh_rcm = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,  // In
                                    lh_sample_point_weights, lh_rc_clipped // In
  );
  stats_update("lh_rcm", lh_rcm, *stats);
  Real1D xtmp    = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "xtmp", ngrdcol);
  auto xtmp_host = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, xtmp);
  for (int i = 0; i < ngrdcol; ++i)
    xtmp_host(i) = vertical_integral(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                     Kokkos::subview(lh_rcm, i, Kokkos::ALL),
                                     Kokkos::subview(gr.dzt, i, Kokkos::ALL));
  Kokkos::deep_copy(xtmp, xtmp_host);
  stats_update("lh_lwp", xtmp, *stats);
  Kokkos::parallel_for(
      "stats_accumulate_lh::weights_sum", ngrdcol, KOKKOS_LAMBDA(int i) {
        real_t value = zero;
        for (int k = 0; k < nzt; ++k)
          for (int sample = 0; sample < num_samples; ++sample)
            value = value + lh_sample_point_weights(i, sample, k);
        xtmp(i) = value;
      });
  stats_update("lh_sample_weights_sum", xtmp, *stats);
  Kokkos::parallel_for(
      "stats_accumulate_lh::weights_average", ngrdcol,
      KOKKOS_LAMBDA(int i) { xtmp(i) = xtmp(i) / real_t((num_samples * nzt)); });
  stats_update("lh_sample_weights_avg", xtmp, *stats);
  const auto lh_thlm = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,   // In
                                           lh_sample_point_weights, lh_thl_clipped // In
  );
  stats_update("lh_thlm", lh_thlm, *stats);
  const auto lh_rvm = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,  // In
                                          lh_sample_point_weights, lh_rv_clipped // In
  );
  stats_update("lh_rvm", lh_rvm, *stats);
  for (int i = 0; i < ngrdcol; ++i)
    xtmp_host(i) = vertical_integral(workspace, nzt, Kokkos::subview(rho_ds_zt, i, Kokkos::ALL),
                                     Kokkos::subview(lh_rvm, i, Kokkos::ALL),
                                     Kokkos::subview(gr.dzt, i, Kokkos::ALL));
  Kokkos::deep_copy(xtmp, xtmp_host);
  stats_update("lh_vwp", xtmp, *stats);
  const auto lh_wm = compute_sample_mean(
      workspace, nzt, num_samples, ngrdcol, // In
      lh_sample_point_weights,
      Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiPDF_w) // In
  );
  stats_update("lh_wm", lh_wm, *stats);
  Real3D lh_hydromet = workspace_scope.view<Real3D>("lh_hydromet", ngrdcol, nzt, hydromet_dim);
  Real4D hydromet_all_points = workspace_scope.view<Real4D>(
      Kokkos::WithoutInitializing, "hydromet_all_points", ngrdcol, num_samples, nzt, hydromet_dim);
  Real3D Ncn_all_points = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing,
                                                       "Ncn_all_points", ngrdcol, num_samples, nzt);
  copy_X_nl_into_hydromet_all_pts(nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs, // In
                                  hydromet_dim, hm_metadata, lh_hydromet,            // In
                                  hydromet_all_points, Ncn_all_points                // Out
  );
  for (int ivar = 0; ivar < hydromet_dim; ++ivar)
    compute_sample_mean(
        workspace, nzt, num_samples, ngrdcol, lh_sample_point_weights,
        Kokkos::subview(hydromet_all_points, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, ivar),
        Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, ivar));
  // Native computes this even if iiPDF_Ncn is absent, in which case its source
  // output was unwritten. Retain that definedness limit rather than zero fill.
  const auto lh_Ncnm = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,   // In
                                           lh_sample_point_weights, Ncn_all_points // In
  );
  stats_update("lh_Ncnm", lh_Ncnm, *stats);
  const auto lh_Ncm = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,  // In
                                          lh_sample_point_weights, lh_Nc_clipped // In
  );
  stats_update("lh_Ncm", lh_Ncm, *stats);
  Real2D lh_cloud_frac =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "lh_cloud_frac", ngrdcol, nzt);
  Kokkos::parallel_for(
      "stats_accumulate_lh::cloud_fraction", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), k = (cell / ngrdcol);
        lh_cloud_frac(i, k) = zero;
        for (int sample = 0; sample < num_samples; ++sample)
          if (X_nl_all_levs(i, sample, k, iiPDF_chi) > zero)
            lh_cloud_frac(i, k) = lh_cloud_frac(i, k) + lh_sample_point_weights(i, sample, k);
        lh_cloud_frac(i, k) = lh_cloud_frac(i, k) / real_t(num_samples);
      });
  stats_update("lh_cloud_frac", lh_cloud_frac, *stats);
  Kokkos::parallel_for(
      "stats_accumulate_lh::cloud_fraction_unweighted", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), k = (cell / ngrdcol);
        lh_cloud_frac(i, k) = zero;
        for (int sample = 0; sample < num_samples; ++sample)
          if (X_nl_all_levs(i, sample, k, iiPDF_chi) > zero)
            lh_cloud_frac(i, k) = lh_cloud_frac(i, k) + one;
        lh_cloud_frac(i, k) = lh_cloud_frac(i, k) / real_t(num_samples);
      });
  stats_update("lh_cloud_frac_unweighted", lh_cloud_frac, *stats);
  const auto lh_chi = compute_sample_mean(
      workspace, nzt, num_samples, ngrdcol, // In
      lh_sample_point_weights,
      Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiPDF_chi) // In
  );
  stats_update("lh_chi", lh_chi, *stats);
  const auto lh_chip2 = compute_sample_variance(
      workspace, nzt, num_samples, ngrdcol, // In
      Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiPDF_chi),
      lh_sample_point_weights, lh_chi // In
  );
  stats_update("lh_chip2", lh_chip2, *stats);
  const auto lh_eta = compute_sample_mean(
      workspace, nzt, num_samples, ngrdcol, // In
      lh_sample_point_weights,
      Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiPDF_eta) // In
  );
  stats_update("lh_eta", lh_eta, *stats);
  const auto lh_wp2_zt = compute_sample_variance(
      workspace, nzt, num_samples, ngrdcol, // In
      Kokkos::subview(X_nl_all_levs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiPDF_w),
      lh_sample_point_weights, lh_wm // In
  );
  stats_update("lh_wp2_zt", lh_wp2_zt, *stats);
  const auto lh_rcp2_zt =
      compute_sample_variance(workspace, nzt, num_samples, ngrdcol,          // In
                              lh_rc_clipped, lh_sample_point_weights, lh_rcm // In
      );
  stats_update("lh_rcp2_zt", lh_rcp2_zt, *stats);
  Real2D total_mean =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "lh_rvm_plus_rcm", ngrdcol, nzt);
  Kokkos::parallel_for(
      "stats_accumulate_lh::total_mean", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ngrdcol), k = (cell / ngrdcol);
        total_mean(i, k) = lh_rvm(i, k) + lh_rcm(i, k);
      });
  const auto lh_rtp2_zt =
      compute_sample_variance(workspace, nzt, num_samples, ngrdcol,              // In
                              lh_rt_clipped, lh_sample_point_weights, total_mean // In
      );
  stats_update("lh_rtp2_zt", lh_rtp2_zt, *stats);
  const auto lh_thlp2_zt =
      compute_sample_variance(workspace, nzt, num_samples, ngrdcol,            // In
                              lh_thl_clipped, lh_sample_point_weights, lh_thlm // In
      );
  stats_update("lh_thlp2_zt", lh_thlp2_zt, *stats);
  if (iirr > -1) {
    const auto lh_rrp2_zt = compute_sample_variance(
        workspace, nzt, num_samples, ngrdcol, // In
        Kokkos::subview(hydromet_all_points, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iirr),
        lh_sample_point_weights, Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iirr) // In
    );
    stats_update("lh_rrp2_zt", lh_rrp2_zt, *stats);
  }
  if (iiPDF_Ncn > -1) {
    const auto lh_Ncnp2_zt =
        compute_sample_variance(workspace, nzt, num_samples, ngrdcol,            // In
                                Ncn_all_points, lh_sample_point_weights, lh_Ncnm // In
        );
    stats_update("lh_Ncnp2_zt", lh_Ncnp2_zt, *stats);
  }
  const auto lh_Ncp2_zt =
      compute_sample_variance(workspace, nzt, num_samples, ngrdcol,          // In
                              lh_Nc_clipped, lh_sample_point_weights, lh_Ncm // In
      );
  stats_update("lh_Ncp2_zt", lh_Ncp2_zt, *stats);
  if (iiNr > -1) {
    const auto lh_Nrp2_zt = compute_sample_variance(
        workspace, nzt, num_samples, ngrdcol, // In
        Kokkos::subview(hydromet_all_points, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iiNr),
        lh_sample_point_weights, Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiNr) // In
    );
    stats_update("lh_Nrp2_zt", lh_Nrp2_zt, *stats);
  }
  // Diagnostic means retain exact Native hydrometeor and number-species order.
  if (iirr > -1)
    stats_update("lh_rrm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iirr), *stats);
  if (iiNr > -1)
    stats_update("lh_Nrm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiNr), *stats);
  if (iiri > -1)
    stats_update("lh_rim", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiri), *stats);
  if (iiNi > -1)
    stats_update("lh_Nim", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiNi), *stats);
  if (iirs > -1)
    stats_update("lh_rsm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iirs), *stats);
  if (iiNs > -1)
    stats_update("lh_Nsm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiNs), *stats);
  if (iirg > -1)
    stats_update("lh_rgm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iirg), *stats);
  if (iiNg > -1)
    stats_update("lh_Ngm", Kokkos::subview(lh_hydromet, Kokkos::ALL, Kokkos::ALL, iiNg), *stats);
}

namespace {

//-----------------------------------------------------------------------------
void stats_accumulate_uniform_lh(
    ScratchWorkspace& workspace,                                                       // InOut
    int nzt, int num_samples, int ngrdcol,                                             // In
    Kokkos::View<std::uint8_t***, Kokkos::LayoutLeft> l_in_precip_all_levs,            // In
    Kokkos::View<int***, Kokkos::LayoutLeft> X_mixt_comp_all_levs,                     // In
    Real3D X_u_chi_all_levs, const pdf_parameter& pdf_params,                          // In
    Real3D lh_sample_point_weights, Kokkos::View<int*, Kokkos::LayoutLeft> k_lh_start, // In
    stats_type& stats                                                                  // InOut
) {
  // Description:
  //   Samples statistics that cannot be deduced from the normal-lognormal
  //   SILHS sample (X_nl_all_levs)
  // References:
  //   None
  // nzt: Number of vertical levels
  // num_samples: Number of SILHS sample points
  // ngrdcol: Number of grid columns
  // l_in_precip_all_levs: Boolean variables indicating whether a sample is in precipitation at a
  //   given height level
  // X_mixt_comp_all_levs: Integers indicating which mixture component a sample is in at a given
  //   height level
  // X_u_chi_all_levs: Uniform value of chi
  // pdf_params: The official PDF parameters!
  // lh_sample_point_weights: The weight of each sample
  // k_lh_start: Zero-based slot for sampling preferentially within       [-] cloud
  auto workspace_scope = workspace.scope();
  if (stats.l_sample) {
    Real3D int_in_precip = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing,
                                                        "int_in_precip", ngrdcol, num_samples, nzt);
    Real3D int_mixt_comp = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing,
                                                        "int_mixt_comp", ngrdcol, num_samples, nzt);
    Real3D one_weights   = workspace_scope.view<Real3D>("one_weights", ngrdcol, num_samples, nzt);
    Kokkos::parallel_for(
        "stats_uniform_lh::precip_indicator", (ngrdcol * num_samples * nzt),
        KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                    k = (cell / (ngrdcol * num_samples));
          if (l_in_precip_all_levs(i, sample, k))
            int_in_precip(i, sample, k) = one;
          else
            int_in_precip(i, sample, k) = zero;
        });
    auto lh_precip_frac = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,  // In
                                              lh_sample_point_weights, int_in_precip // In
    );
    stats_update("lh_precip_frac", lh_precip_frac, stats);
    Kokkos::deep_copy(one_weights, one);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, one_weights, int_in_precip,
                        lh_precip_frac);
    stats_update("lh_precip_frac_unweighted", lh_precip_frac, stats);
    Kokkos::parallel_for(
        "stats_uniform_lh::mixture_indicator", (ngrdcol * num_samples * nzt),
        KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), sample = ((cell / ngrdcol) % num_samples),
                    k = (cell / (ngrdcol * num_samples));
          if (X_mixt_comp_all_levs(i, sample, k) == 1)
            int_mixt_comp(i, sample, k) = one;
          else
            int_mixt_comp(i, sample, k) = zero;
        });
    auto lh_mixt_frac = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,  // In
                                            lh_sample_point_weights, int_mixt_comp // In
    );
    stats_update("lh_mixt_frac", lh_mixt_frac, stats);
    compute_sample_mean(workspace, nzt, num_samples, ngrdcol, one_weights, int_mixt_comp,
                        lh_mixt_frac);
    stats_update("lh_mixt_frac_unweighted", lh_mixt_frac, stats);
    Real1D k_start_real =
        workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "k_lh_start_real", ngrdcol);
    Kokkos::parallel_for(
        "stats_uniform_lh::height_index", ngrdcol,
        KOKKOS_LAMBDA(int i) { k_start_real(i) = real_t(k_lh_start(i) + 1); });
    stats_update("k_lh_start", k_start_real, stats);
    const auto importance_categories = define_importance_categories();
    const auto cloud_frac_1 = pdf_params.cloud_frac_1, cloud_frac_2 = pdf_params.cloud_frac_2;
    Real3D lh_samp_frac = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing, "lh_samp_frac",
                                                       ngrdcol, nzt, num_importance_categories);
    Kokkos::parallel_for(
        "stats_uniform_lh::categories", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
          const int i = (cell % ngrdcol), k = (cell / ngrdcol);
          Kokkos::Array<int, num_importance_categories> category_counts;
          for (int j = 0; j < num_importance_categories; ++j)
            category_counts[j] = 0;
          for (int isample = 0; isample < num_samples; ++isample) {
            bool l_in_comp_1;
            real_t cloud_frac_i;
            if (X_mixt_comp_all_levs(i, isample, k) == 1) {
              l_in_comp_1  = true;
              cloud_frac_i = cloud_frac_1(i, k);
            } else {
              l_in_comp_1  = false;
              cloud_frac_i = cloud_frac_2(i, k);
            }
            const bool l_in_cloud = X_u_chi_all_levs(i, isample, k) > one - cloud_frac_i;
            for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
              if (l_in_cloud == importance_categories[icategory].l_in_cloud &&
                  bool(l_in_precip_all_levs(i, isample, k)) ==
                      importance_categories[icategory].l_in_precip &&
                  l_in_comp_1 == importance_categories[icategory].l_in_component_1) {
                ++category_counts[icategory];
                break;
              }
            }
          }
          for (int icategory = 0; icategory < num_importance_categories; ++icategory)
            lh_samp_frac(i, k, icategory) =
                real_t(category_counts[icategory]) / real_t(num_samples);
          // Native suppresses array index1 even on a descending grid.
          if (k == 0)
            for (int icategory = 0; icategory < num_importance_categories; ++icategory)
              lh_samp_frac(i, k, icategory) = zero;
        });
    for (int icategory = 0; icategory < num_importance_categories; ++icategory) {
      const std::string samp_frac_name = "lh_samp_frac_" + std::to_string((icategory + 1));
      stats_update(samp_frac_name,
                   Kokkos::subview(lh_samp_frac, Kokkos::ALL, Kokkos::ALL, icategory), stats);
    }
  }
}
} // namespace

//-----------------------------------------------------------------------------
void copy_X_nl_into_hydromet_all_pts(int nzt, int pdf_dim, int num_samples, int ngrdcol,
                                     Real4D X_nl_all_levs, // In
                                     int hydromet_dim, const hm_metadata_type& hm_metadata,
                                     Real3D hydromet,                                  // In
                                     Real4D hydromet_all_points, Real3D Ncn_all_points // Out
) {
  // Copy each sampled hydrometeor into its species array; unsampled species
  // retain the mean field. Metadata indices and mixture IDs remain one based.
  // hydromet_dim: Number of hydrometeor species
  // hydromet: Hydrometeor species    [units vary]
  for (int sample = 0; sample < num_samples; ++sample) {
    for (int ivar = 0; ivar < hydromet_dim; ++ivar) {
      int pdf_hydromet_idx = -2;
      if (ivar == hm_metadata.iirr && hm_metadata.iiPDF_rr > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_rr;
      else if (ivar == hm_metadata.iirs && hm_metadata.iiPDF_rs > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_rs;
      else if (ivar == hm_metadata.iiri && hm_metadata.iiPDF_ri > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_ri;
      else if (ivar == hm_metadata.iirg && hm_metadata.iiPDF_rg > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_rg;
      else if (ivar == hm_metadata.iiNr && hm_metadata.iiPDF_Nr > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_Nr;
      else if (ivar == hm_metadata.iiNs && hm_metadata.iiPDF_Ns > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_Ns;
      else if (ivar == hm_metadata.iiNg && hm_metadata.iiPDF_Ng > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_Ng;
      else if (ivar == hm_metadata.iiNi && hm_metadata.iiPDF_Ni > -1)
        pdf_hydromet_idx = hm_metadata.iiPDF_Ni;
      if (pdf_hydromet_idx > -1) {
        Kokkos::parallel_for(
            "copy_X_nl::sampled", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
              const int i = (cell % ngrdcol);
              const int k = (cell / ngrdcol);
              hydromet_all_points(i, sample, k, ivar) =
                  real_t(X_nl_all_levs(i, sample, k, pdf_hydromet_idx));
            });
      } else {
        Kokkos::parallel_for(
            "copy_X_nl::mean", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
              const int i                             = (cell % ngrdcol);
              const int k                             = (cell / ngrdcol);
              hydromet_all_points(i, sample, k, ivar) = hydromet(i, k, ivar);
            });
      }
    }
    // Native leaves Ncn output unwritten when the variate is absent. Do not
    // silently adopt JAX's zero fill or assume defined downstream diagnostics.
    const int iiPDF_Ncn = hm_metadata.iiPDF_Ncn;
    if (iiPDF_Ncn > -1) {
      Kokkos::parallel_for(
          "copy_X_nl::Ncn", (ngrdcol * nzt), KOKKOS_LAMBDA(int cell) {
            const int i                  = (cell % ngrdcol);
            const int k                  = (cell / ngrdcol);
            Ncn_all_points(i, sample, k) = real_t(X_nl_all_levs(i, sample, k, iiPDF_Ncn));
          });
    }
  }
}
} // namespace clubb_kokkos
