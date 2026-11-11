// SILHS microphysical source terms for variances and covariances.
// Public function: lh_microphys_var_covar_driver_api.
#include "scratch_workspace.hpp"
#include "SILHS/lh_microphys_var_covar_module.hpp"
#include "SILHS/math_utilities.hpp"
#include "CLUBB_core/constants_clubb.hpp"

namespace clubb_kokkos {

//-----------------------------------------------------------------------
void lh_microphys_var_covar_driver_api(
    ScratchWorkspace& workspace,                                                      // InOut
    int nzt, int num_samples, int ngrdcol, real_t dt, Real3D lh_sample_point_weights, // In
    const pdf_parameter& pdf_params, Real3D lh_rt_all, Real3D lh_thl_all,             // In
    Real3D lh_w_all,                                                                  // In
    Real3D lh_rcm_mc_all, Real3D lh_rvm_mc_all, Real3D lh_thlm_mc_all,                // In
    bool l_lh_instant_var_covar_src,                                                  // In
    Real2D lh_rtp2_mc_zt, Real2D lh_thlp2_mc_zt, Real2D lh_wprtp_mc_zt,               // Out
    Real2D lh_wpthlp_mc_zt, Real2D lh_rtpthlp_mc_zt                                   // Out
) {
  // Compute microphysical tendencies in all columns. Explicit time integration
  // gives rtp2_mc=2*cov(rt,rt_mc)+dt*var(rt_mc), similarly for thlp2; the
  // covariance tendency adds dt*cov(rt_mc,thl_mc). The instantaneous option
  // omits these dt terms to agree with the KK upscaled analytic solution.
  // See clubb:ticket:753.
  // nzt: Number of vertical levels
  // num_samples: Number of SILHS sample points
  // ngrdcol: Number of model columns
  // dt: Model time step                             [s]
  // lh_sample_point_weights: Weight of SILHS sample points
  // pdf_params: The PDF parameters_silhs
  // lh_rt_all: SILHS samples of total water                [kg/kg]
  // lh_thl_all: SILHS samples of potential temperature      [K]
  // lh_w_all: SILHS samples of vertical velocity          [m/s]
  // lh_rcm_mc_all: SILHS microphys. tendency of rcm            [kg/kg/s]
  // lh_rvm_mc_all: SILHS microphys. tendency of rvm            [kg/kg/s]
  // lh_thlm_mc_all: SILHS microphys. tendency of thlm           [K/s]
  // l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies  [-]
  // lh_rtp2_mc_zt: SILHS microphys. est. tendency of <rt'^2>   [(kg/kg)^2/s]
  // lh_thlp2_mc_zt: SILHS microphys. est. tendency of <thl'^2>  [K^2/s]
  // lh_wprtp_mc_zt: SILHS microphys. est. tendency of <w'rt'>   [m*(kg/kg)/s^2]
  // lh_wpthlp_mc_zt: SILHS microphys. est. tendency of <w'thl'>  [m*K/s^2]
  // lh_rtpthlp_mc_zt: SILHS microphys. est. tendency of <rt'thl'> [K*(kg/kg)/s]
  auto workspace_scope = workspace.scope();
  const auto mixt_frac = pdf_params.mixt_frac;
  const auto rt_1      = pdf_params.rt_1;
  const auto rt_2      = pdf_params.rt_2;
  const auto thl_1     = pdf_params.thl_1;
  const auto thl_2     = pdf_params.thl_2;
  const auto w_1       = pdf_params.w_1;
  const auto w_2       = pdf_params.w_2;
  Real3D lh_rt_mc_all  = workspace_scope.view<Real3D>(Kokkos::WithoutInitializing, "lh_rt_mc_all",
                                                      ngrdcol, num_samples, nzt);
  Real2D mean_rt =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "mean_rt", ngrdcol, nzt);
  Real2D mean_thl =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "mean_thl", ngrdcol, nzt);
  Real2D mean_w = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "mean_w", ngrdcol, nzt);
  Kokkos::parallel_for(
      "lh_microphys_var_covar::rt_mc", ngrdcol * num_samples * nzt, KOKKOS_LAMBDA(int cell) {
        const int i                = cell % ngrdcol;
        const int sample           = (cell / ngrdcol) % num_samples;
        const int k                = cell / (ngrdcol * num_samples);
        lh_rt_mc_all(i, sample, k) = lh_rcm_mc_all(i, sample, k) + lh_rvm_mc_all(i, sample, k);
      });
  // Calculate means, variances and covariances needed for the tendency terms.
  Kokkos::parallel_for(
      "lh_microphys_var_covar::means", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i    = cell % ngrdcol;
        const int k    = cell / ngrdcol;
        mean_rt(i, k)  = mixt_frac(i, k) * rt_1(i, k) + (one - mixt_frac(i, k)) * rt_2(i, k);
        mean_thl(i, k) = mixt_frac(i, k) * thl_1(i, k) + (one - mixt_frac(i, k)) * thl_2(i, k);
        mean_w(i, k)   = mixt_frac(i, k) * w_1(i, k) + (one - mixt_frac(i, k)) * w_2(i, k);
      });
  const auto mean_rt_mc = compute_sample_mean(workspace, nzt, num_samples, ngrdcol, // In
                                              lh_sample_point_weights, lh_rt_mc_all // In
  );
  const auto covar_rt_rt_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,        // In
                                lh_sample_point_weights, lh_rt_all, mean_rt, // In
                                lh_rt_mc_all, mean_rt_mc                     // In
      );
  const auto mean_thl_mc = compute_sample_mean(workspace, nzt, num_samples, ngrdcol,   // In
                                               lh_sample_point_weights, lh_thlm_mc_all // In
  );
  const auto covar_thl_thl_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,          // In
                                lh_sample_point_weights, lh_thl_all, mean_thl, // In
                                lh_thlm_mc_all, mean_thl_mc                    // In
      );
  const auto covar_w_rt_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,      // In
                                lh_sample_point_weights, lh_w_all, mean_w, // In
                                lh_rt_mc_all, mean_rt_mc                   // In
      );
  const auto covar_w_thl_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,      // In
                                lh_sample_point_weights, lh_w_all, mean_w, // In
                                lh_thlm_mc_all, mean_thl_mc                // In
      );
  const auto covar_thl_rt_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,          // In
                                lh_sample_point_weights, lh_thl_all, mean_thl, // In
                                lh_rt_mc_all, mean_rt_mc                       // In
      );
  const auto covar_rt_thl_mc =
      compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,        // In
                                lh_sample_point_weights, lh_rt_all, mean_rt, // In
                                lh_thlm_mc_all, mean_thl_mc                  // In
      );
  // Compute optional corrections before the output kernel so each tendency is
  // completed in one traversal, preserving the base-plus-dt operation order.
  Real2D var_rt_mc, var_thl_mc, covar_rt_mc_thl_mc;
  if (!l_lh_instant_var_covar_src) {
    // These timestep-dependent terms arise from explicitly integrated means;
    // they are absent from KK upscaled and cause non-convergence with it.
    var_rt_mc  = compute_sample_variance(workspace, nzt, num_samples, ngrdcol,             // In
                                         lh_rt_mc_all, lh_sample_point_weights, mean_rt_mc // In
    );
    var_thl_mc = compute_sample_variance(workspace, nzt, num_samples, ngrdcol,                // In
                                         lh_thlm_mc_all, lh_sample_point_weights, mean_thl_mc // In
    );
    covar_rt_mc_thl_mc =
        compute_sample_covariance(workspace, nzt, num_samples, ngrdcol,              // In
                                  lh_sample_point_weights, lh_rt_mc_all, mean_rt_mc, // In
                                  lh_thlm_mc_all, mean_thl_mc                        // In
        );
  }
  Kokkos::parallel_for(
      "lh_microphys_var_covar::tendencies", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i            = cell % ngrdcol;
        const int k            = cell / ngrdcol;
        lh_rtp2_mc_zt(i, k)    = two * covar_rt_rt_mc(i, k);
        lh_thlp2_mc_zt(i, k)   = two * covar_thl_thl_mc(i, k);
        lh_wprtp_mc_zt(i, k)   = covar_w_rt_mc(i, k);
        lh_wpthlp_mc_zt(i, k)  = covar_w_thl_mc(i, k);
        lh_rtpthlp_mc_zt(i, k) = covar_thl_rt_mc(i, k) + covar_rt_thl_mc(i, k);
        if (!l_lh_instant_var_covar_src) {
          lh_rtp2_mc_zt(i, k)    = lh_rtp2_mc_zt(i, k) + dt * var_rt_mc(i, k);
          lh_thlp2_mc_zt(i, k)   = lh_thlp2_mc_zt(i, k) + dt * var_thl_mc(i, k);
          lh_rtpthlp_mc_zt(i, k) = lh_rtpthlp_mc_zt(i, k) + dt * covar_rt_mc_thl_mc(i, k);
        }
      });
}

} // namespace clubb_kokkos
