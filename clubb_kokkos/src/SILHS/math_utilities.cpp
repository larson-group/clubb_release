// Weighted SILHS sample moments and the Native signed MT95 integer draw.
// Public functions: compute_sample_mean, compute_sample_variance,
// compute_sample_covariance and rand_integer_in_range.
#include "SILHS/math_utilities.hpp"
#include "SILHS/mt95.hpp"
#include <cstdlib>

namespace clubb_kokkos {

template <class Source> static bool overlapping_fields(const Source& source, Real2D output) {
  // Compare complete same-space spans so aliases retain input snapshot semantics.
  if (source.size() == 0 || output.size() == 0)
    return false;
  const auto source_begin = reinterpret_cast<std::uintptr_t>(source.data());
  const auto output_begin = reinterpret_cast<std::uintptr_t>(output.data());
  if (source_begin <= output_begin)
    return output_begin - source_begin < source.span() * sizeof(real_t);
  return source_begin - output_begin < output.span() * sizeof(real_t);
}

static void compute_sample_mean_values(int n_levels, int n_samples, int ngrdcol, Real3D weight,
                                       Real3D x_sample, Real2D mean) {
  // Find the weighted mean of sample points in every model column.
  // Each independent cell retains Native's increasing sample accumulation.
  // n_levels: Number of sample levels
  // n_samples: Number of sample points
  // ngrdcol: Number of model columns
  // weight: Weights for individual points of the vector
  // x_sample: Collection of sample points    [units vary]
  // Keep the increasing sample fold in a scalar; store the final field once.
  Kokkos::parallel_for(
      "compute_sample_mean", ngrdcol * n_levels, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t mean_sum = real_t(0);
        for (int sample = 0; sample < n_samples; ++sample)
          mean_sum = mean_sum + weight(i, sample, k) * x_sample(i, sample, k);
        mean(i, k) = mean_sum / real_t(n_samples);
      });
}

static void compute_sample_variance_values(int n_levels, int n_samples, int ngrdcol,
                                           Real3D x_sample, Real3D weight, Real2D x_mean,
                                           Real2D variance) {
  // Compute weighted population variance; normalization is N, not N-1.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the variance of
  // ngrdcol: Number of model columns
  // x_sample: Collection of sample points    [units vary]
  // weight: Coefficient to weight the nth sample point by [-]
  // x_mean: Mean sample points [units vary]
  // Keep the increasing sample fold in a scalar; store the final field once.
  Kokkos::parallel_for(
      "compute_sample_variance", ngrdcol * n_levels, KOKKOS_LAMBDA(int cell) {
        const int i    = cell % ngrdcol;
        const int k    = cell / ngrdcol;
        real_t variance_sum = real_t(0);
        for (int sample = 0; sample < n_samples; ++sample) {
          const real_t x_deviation = x_sample(i, sample, k) - x_mean(i, k);
          variance_sum = variance_sum + weight(i, sample, k) * (x_deviation * x_deviation);
        }
        variance(i, k) = variance_sum / real_t(n_samples);
      });
}

static void compute_sample_covariance_values(int n_levels, int n_samples, int ngrdcol,
                                             Real3D weight, Real3D x_sample, Real2D x_mean,
                                             Real3D y_sample, Real2D y_mean, Real2D covariance) {
  // Compute weighted covariance of two sample collections in every column.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the covariance of
  // ngrdcol: Number of model columns
  // weight: Coefficient to weight the nth sample point by [-]
  // x_sample: Collection of sample points    [units vary]
  // x_mean: Mean sample points [units vary]
  // y_sample: Collection of sample points    [units vary]
  // y_mean: Mean sample points [units vary]
  // Keep the increasing sample fold in a scalar; store the final field once.
  Kokkos::parallel_for(
      "compute_sample_covariance", ngrdcol * n_levels, KOKKOS_LAMBDA(int cell) {
        const int i      = cell % ngrdcol;
        const int k      = cell / ngrdcol;
        real_t covariance_sum = real_t(0);
        for (int sample = 0; sample < n_samples; ++sample) {
          const real_t x_deviation = x_sample(i, sample, k) - x_mean(i, k);
          const real_t y_deviation = y_sample(i, sample, k) - y_mean(i, k);
          covariance_sum = covariance_sum + weight(i, sample, k) * x_deviation * y_deviation;
        }
        covariance(i, k) = covariance_sum / real_t(n_samples);
      });
}

Real2D compute_sample_mean(int n_levels, int n_samples, int ngrdcol, Real3D weight,
                           Real3D x_sample) {
  // Find the weighted mean of sample points in every model column.
  // Each independent cell retains Native's increasing sample accumulation.
  // n_levels: Number of sample levels
  // n_samples: Number of sample points
  // ngrdcol: Number of model columns
  // weight: Weights for individual points of the vector
  // x_sample: Collection of sample points    [units vary]
  Real2D mean(Kokkos::view_alloc(Kokkos::WithoutInitializing, "mean"), ngrdcol, n_levels);
  compute_sample_mean_values(n_levels, n_samples, ngrdcol, weight, x_sample, mean);
  return mean;
}

Real2D compute_sample_variance(int n_levels, int n_samples, int ngrdcol, Real3D x_sample,
                               Real3D weight, Real2D x_mean) {
  // Compute weighted population variance; normalization is N, not N-1.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the variance of
  // ngrdcol: Number of model columns
  // x_sample: Collection of sample points    [units vary]
  // weight: Coefficient to weight the nth sample point by [-]
  // x_mean: Mean sample points [units vary]
  Real2D variance(Kokkos::view_alloc(Kokkos::WithoutInitializing, "variance"), ngrdcol, n_levels);
  compute_sample_variance_values(n_levels, n_samples, ngrdcol, x_sample, weight, x_mean, variance);
  return variance;
}

Real2D compute_sample_covariance(int n_levels, int n_samples, int ngrdcol, Real3D weight,
                                 Real3D x_sample, Real2D x_mean, Real3D y_sample, Real2D y_mean) {
  // Compute weighted covariance of two sample collections in every column.
  // n_levels: Number of sample levels in the mean / variance
  // n_samples: Number of sample points to compute the covariance of
  // ngrdcol: Number of model columns
  // weight: Coefficient to weight the nth sample point by [-]
  // x_sample: Collection of sample points    [units vary]
  // x_mean: Mean sample points [units vary]
  // y_sample: Collection of sample points    [units vary]
  // y_mean: Mean sample points [units vary]
  Real2D covariance(Kokkos::view_alloc(Kokkos::WithoutInitializing, "covariance"), ngrdcol,
                    n_levels);
  compute_sample_covariance_values(n_levels, n_samples, ngrdcol, weight, x_sample, x_mean, y_sample,
                                   y_mean, covariance);
  return covariance;
}

Real2D compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                           Real3D weight, Real3D x_sample) {
  // Return a fully written field reserved in the active caller frame.
  auto mean = workspace.view<Real2D>(Kokkos::WithoutInitializing, "mean", ngrdcol, n_levels);
  compute_sample_mean_values(n_levels, n_samples, ngrdcol, weight, x_sample, mean);
  return mean;
}

void compute_sample_mean(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                         Real3D weight, Real3D x_sample, Real2D mean) {
  // Preserve the owned-return snapshot when the explicit output overlaps an input.
  if (overlapping_fields(weight, mean) || overlapping_fields(x_sample, mean)) {
    auto overlap_scope = workspace.scope();
    auto result = overlap_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                             "compute_sample_mean::overlap", ngrdcol, n_levels);
    compute_sample_mean_values(n_levels, n_samples, ngrdcol, weight, x_sample, result);
    Kokkos::deep_copy(mean, result);
  } else {
    compute_sample_mean_values(n_levels, n_samples, ngrdcol, weight, x_sample, mean);
  }
}

Real2D compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Real3D x_sample, Real3D weight, Real2D x_mean) {
  // Return a fully written field reserved in the active caller frame.
  auto variance =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "variance", ngrdcol, n_levels);
  compute_sample_variance_values(n_levels, n_samples, ngrdcol, x_sample, weight, x_mean, variance);
  return variance;
}

void compute_sample_variance(ScratchWorkspace& workspace, int n_levels, int n_samples, int ngrdcol,
                             Real3D x_sample, Real3D weight, Real2D x_mean, Real2D variance) {
  // Preserve the owned-return snapshot when the explicit output overlaps an input.
  if (overlapping_fields(x_sample, variance) || overlapping_fields(weight, variance) ||
      overlapping_fields(x_mean, variance)) {
    auto overlap_scope = workspace.scope();
    auto result = overlap_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                             "compute_sample_variance::overlap", ngrdcol, n_levels);
    compute_sample_variance_values(n_levels, n_samples, ngrdcol, x_sample, weight, x_mean, result);
    Kokkos::deep_copy(variance, result);
  } else {
    compute_sample_variance_values(n_levels, n_samples, ngrdcol, x_sample, weight, x_mean,
                                   variance);
  }
}

Real2D compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                                 int ngrdcol, Real3D weight, Real3D x_sample, Real2D x_mean,
                                 Real3D y_sample, Real2D y_mean) {
  // Return a fully written field reserved in the active caller frame.
  auto covariance =
      workspace.view<Real2D>(Kokkos::WithoutInitializing, "covariance", ngrdcol, n_levels);
  compute_sample_covariance_values(n_levels, n_samples, ngrdcol, weight, x_sample, x_mean, y_sample,
                                   y_mean, covariance);
  return covariance;
}

void compute_sample_covariance(ScratchWorkspace& workspace, int n_levels, int n_samples,
                               int ngrdcol, Real3D weight, Real3D x_sample, Real2D x_mean,
                               Real3D y_sample, Real2D y_mean, Real2D covariance) {
  // Preserve the owned-return snapshot when the explicit output overlaps an input.
  if (overlapping_fields(weight, covariance) || overlapping_fields(x_sample, covariance) ||
      overlapping_fields(x_mean, covariance) || overlapping_fields(y_sample, covariance) ||
      overlapping_fields(y_mean, covariance)) {
    auto overlap_scope = workspace.scope();
    auto result        = overlap_scope.view<Real2D>(
        Kokkos::WithoutInitializing, "compute_sample_covariance::overlap", ngrdcol, n_levels);
    compute_sample_covariance_values(n_levels, n_samples, ngrdcol, weight, x_sample, x_mean,
                                     y_sample, y_mean, result);
    Kokkos::deep_copy(covariance, result);
  } else {
    compute_sample_covariance_values(n_levels, n_samples, ngrdcol, weight, x_sample, x_mean,
                                     y_sample, y_mean, covariance);
  }
}

//-----------------------------------------------------------------------
int rand_integer_in_range(int low, int high, // In
                          mt95& rng          // Case state
) {
  // Return a signed-MT95 draw in [low,high]. Native MOD introduces a slight
  // bias toward smaller integers; retain that algorithm and its draw count.
  // low: Lowest possible returned value
  // high: Highest possible returned value
  const int range_width = high - low + 1;
  std::int32_t rand_32;
  rng.genrand_int32(rand_32);
  return std::abs(rand_32 % range_width) + low;
}

} // namespace clubb_kokkos
