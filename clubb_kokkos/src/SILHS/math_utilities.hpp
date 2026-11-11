// Weighted SILHS sample means, variances, covariances and signed RNG indices.
// Public functions: compute_sample_mean, compute_sample_variance,
// compute_sample_covariance and rand_integer_in_range. Sample folds stay ordered.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include "SILHS/mt95.hpp"

namespace clubb_kokkos {
int rand_integer_in_range(int low, int high, // In
                          mt95& rng          // Case state
);
Real2D compute_sample_mean(int n_levels, int n_samples, int ngrdcol, // In
                           Real3D weight, Real3D x_sample            // In
);

Real2D compute_sample_variance(int n_levels, int n_samples, int ngrdcol,     // In
                               Real3D x_sample, Real3D weight, Real2D x_mean // In
);

Real2D compute_sample_covariance(int n_levels, int n_samples, int ngrdcol,      // In
                                 Real3D weight, Real3D x_sample, Real2D x_mean, // In
                                 Real3D y_sample, Real2D y_mean                 // In
);

// Scoped results live in the active caller frame. Explicit outputs retain their
// allocation and may overlap inputs; owning interfaces above remain available.
Real2D compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                           Real3D weight, Real3D x_sample);
void compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                         Real3D weight, Real3D x_sample, Real2D mean);
Real2D compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Real3D x_sample, Real3D weight, Real2D x_mean);
void compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                             Real3D x_sample, Real3D weight, Real2D x_mean, Real2D variance);
Real2D compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                                 int ngrdcol, Real3D weight, Real3D x_sample, Real2D x_mean,
                                 Real3D y_sample, Real2D y_mean);
void compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Real3D weight, Real3D x_sample, Real2D x_mean,
                               Real3D y_sample, Real2D y_mean, Real2D covariance);

template <class Field_weight, class Field_x_sample>
Real2D compute_sample_mean(int n_levels, int n_samples, int ngrdcol,    // In
                           Field_weight weight, Field_x_sample x_sample // In
) {
  // Shared-owner conversion at the host/strided public field boundary.
  // n_levels: Number of sample levels
  // n_samples: Number of sample points
  // ngrdcol: Number of model columns
  // weight: Weights for individual points of the vector
  // x_sample: Collection of sample points    [units vary]
  return compute_sample_mean(n_levels, n_samples, ngrdcol,                            // In
                             to_execution_space(weight), to_execution_space(x_sample) // In
  );
}

template <class Field_x_sample, class Field_weight, class Field_x_mean>
Real2D compute_sample_variance(int n_levels, int n_samples, int ngrdcol, // In
                               Field_x_sample x_sample, Field_weight weight,
                               Field_x_mean x_mean // In
) {
  // Shared-owner conversion at the host/strided public field boundary.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the variance of
  // ngrdcol: Number of model columns
  // x_sample: Collection of sample points    [units vary]
  // weight: Coefficient to weight the nth sample point by [-]
  // x_mean: Mean sample points [units vary]
  return compute_sample_variance(n_levels, n_samples, ngrdcol, // In
                                 to_execution_space(x_sample), to_execution_space(weight),
                                 to_execution_space(x_mean) // In
  );
}

template <class Field_weight, class Field_x_sample, class Field_x_mean, class Field_y_sample,
          class Field_y_mean>
Real2D compute_sample_covariance(int n_levels, int n_samples, int ngrdcol, // In
                                 Field_weight weight, Field_x_sample x_sample,
                                 Field_x_mean x_mean,                         // In
                                 Field_y_sample y_sample, Field_y_mean y_mean // In
) {
  // Shared-owner conversion at the host/strided public field boundary.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the covariance of
  // ngrdcol: Number of model columns
  // weight: Coefficient to weight the nth sample point by [-]
  // x_sample: Collection of sample points    [units vary]
  // x_mean: Mean sample points [units vary]
  // y_sample: Collection of sample points    [units vary]
  // y_mean: Mean sample points [units vary]
  return compute_sample_covariance(n_levels, n_samples, ngrdcol, // In
                                   to_execution_space(weight), to_execution_space(x_sample),
                                   to_execution_space(x_mean),                              // In
                                   to_execution_space(y_sample), to_execution_space(y_mean) // In
  );
}

template <class Field_weight, class Field_x_sample>
Real2D compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                           Field_weight weight, Field_x_sample x_sample) {
  // Stage transient inputs in the caller frame; the returned execution field outlives this helper.
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  return compute_sample_mean(workspace, n_levels, n_samples, ngrdcol, weight_execution,
                             x_sample_execution);
}

template <class Field_weight, class Field_x_sample, class Output>
void compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                         Field_weight weight, Field_x_sample x_sample, Output mean) {
  // Fill existing host/strided storage; staging ends after output restoration.
  auto staging_scope      = workspace.scope();
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  auto result             = to_execution_space(workspace, Kokkos::WithoutInitializing, mean);
  compute_sample_mean(workspace, n_levels, n_samples, ngrdcol, weight_execution, x_sample_execution,
                      result);
  from_execution_space(workspace, mean, result);
}

template <class Field_x_sample, class Field_weight, class Field_x_mean>
Real2D compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Field_x_sample x_sample, Field_weight weight,
                               Field_x_mean x_mean) {
  // Stage transient inputs in the caller frame; the returned execution field outlives this helper.
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_mean_execution   = to_execution_space(workspace, x_mean);
  return compute_sample_variance(workspace, n_levels, n_samples, ngrdcol, x_sample_execution,
                                 weight_execution, x_mean_execution);
}

template <class Field_x_sample, class Field_weight, class Field_x_mean, class Output>
void compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                             Field_x_sample x_sample, Field_weight weight, Field_x_mean x_mean,
                             Output variance) {
  // Fill existing host/strided storage; staging ends after output restoration.
  auto staging_scope      = workspace.scope();
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_mean_execution   = to_execution_space(workspace, x_mean);
  auto result             = to_execution_space(workspace, Kokkos::WithoutInitializing, variance);
  compute_sample_variance(workspace, n_levels, n_samples, ngrdcol, x_sample_execution,
                          weight_execution, x_mean_execution, result);
  from_execution_space(workspace, variance, result);
}

template <class Field_weight, class Field_x_sample, class Field_x_mean, class Field_y_sample,
          class Field_y_mean>
Real2D compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                                 int ngrdcol, Field_weight weight, Field_x_sample x_sample,
                                 Field_x_mean x_mean, Field_y_sample y_sample,
                                 Field_y_mean y_mean) {
  // Stage transient inputs in the caller frame; the returned execution field outlives this helper.
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  auto x_mean_execution   = to_execution_space(workspace, x_mean);
  auto y_sample_execution = to_execution_space(workspace, y_sample);
  auto y_mean_execution   = to_execution_space(workspace, y_mean);
  return compute_sample_covariance(workspace, n_levels, n_samples, ngrdcol, weight_execution,
                                   x_sample_execution, x_mean_execution, y_sample_execution,
                                   y_mean_execution);
}

template <class Field_weight, class Field_x_sample, class Field_x_mean, class Field_y_sample,
          class Field_y_mean, class Output>
void compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Field_weight weight, Field_x_sample x_sample,
                               Field_x_mean x_mean, Field_y_sample y_sample, Field_y_mean y_mean,
                               Output covariance) {
  // Fill existing host/strided storage; staging ends after output restoration.
  auto staging_scope      = workspace.scope();
  auto weight_execution   = to_execution_space(workspace, weight);
  auto x_sample_execution = to_execution_space(workspace, x_sample);
  auto x_mean_execution   = to_execution_space(workspace, x_mean);
  auto y_sample_execution = to_execution_space(workspace, y_sample);
  auto y_mean_execution   = to_execution_space(workspace, y_mean);
  auto result             = to_execution_space(workspace, Kokkos::WithoutInitializing, covariance);
  compute_sample_covariance(workspace, n_levels, n_samples, ngrdcol, weight_execution,
                            x_sample_execution, x_mean_execution, y_sample_execution,
                            y_mean_execution, result);
  from_execution_space(workspace, covariance, result);
}

} // namespace clubb_kokkos
