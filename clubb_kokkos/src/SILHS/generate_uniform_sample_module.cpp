// Ordered uniform Latin-hypercube draws and Native Fisher-Yates permutations.
#include "SILHS/generate_uniform_sample_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include <iostream>
#include <limits>

namespace clubb_kokkos {
//-----------------------------------------------------------------------
real_t rand_uniform_real(mt95& rng // Case state
) {
  // Generate a uniform real in (0,1) with the Native signed MT95 conversion.
  double rand_uniform_real_genrand_real;
  rng.genrand_real3(rand_uniform_real_genrand_real);
  real_t value = real_t(rand_uniform_real_genrand_real);
  // Kind conversion can round a Native open-interval draw to one.
  if (value >= one)
    value = one - std::numeric_limits<real_t>::epsilon();
  return value;
}

//-------------------------------------------------------------------------------
void generate_uniform_lh_sample(int iter, int num_samples, int sequence_length, int n_vars, // In
                                bool l_lh_deterministic_test,                               // In
                                Real2D X_u_one_lev,                                         // Out
                                latin_hypercube_arrays_state& sampling_state, // Case state
                                mt95& rng                                     // Case state
) {
  // Generate a uniform Latin-hypercube matrix (sample,variate), preserving the
  // source permutation between refreshes and its first-num_samples row choice.
  // Art B. Owen (2003), Quasi-Monte Carlo Sampling, SIGGRAPH2003;
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:lh_algorithm.
  // iter: Model iteration number
  // num_samples: `n'  Number of samples generated
  // sequence_length: `n_t' Num. random samples before sequence repeats
  // n_vars: Number of uniform variables to generate
  // l_lh_deterministic_test: Ordered strata and repeating draws for repeatable testing
  // X_u_one_lev: num_samples by n_vars matrix, X each row of which is a n_vars-dimensional sample
  constexpr real_t test_uniform_draws[4] = {real_t(.125), real_t(.625), real_t(.375), real_t(.875)};
  constexpr bool l_diagnostic_iter_check = true;
  const int nt_repeat                    = num_samples * sequence_length;
  if (!sampling_state.one_height_time_matrix.data()) {
    sampling_state.one_height_time_matrix =
        SilhsInt2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "one_height_time_matrix"),
                   nt_repeat, n_vars);
    sampling_state.prior_iter = iter;
  } else if (l_diagnostic_iter_check && sequence_length > 1) {
    // Duplicate calls warn without changing prior_iter; a valid next iteration
    // updates it. This also preserves multi-column Native diagnostics.
    if (sampling_state.prior_iter != iter - 1)
      std::cerr << "The iteration number in latin_hypercube_driver is not incrementing properly.\n";
    else
      sampling_state.prior_iter = iter;
  }
  const int i_rmd = (iter - 1) % sequence_length;
  if (i_rmd == 0) {
    if (l_lh_deterministic_test) {
      auto host = Kokkos::create_mirror_view(sampling_state.one_height_time_matrix);
      for (int k = 0; k < n_vars; ++k)
        for (int j = 0; j < nt_repeat; ++j)
          host(j, k) = j;
      Kokkos::deep_copy(sampling_state.one_height_time_matrix, host);
    } else
      permute_height_time(nt_repeat, n_vars, sampling_state.one_height_time_matrix, // In / Out
                          rng                                                       // Case state
      );
  }
  const auto matrix = to_host_space(sampling_state.one_height_time_matrix);
  auto values       = Kokkos::create_mirror_view(X_u_one_lev);
  // Draw order is sample then variate; RNG state belongs to the caller.
  for (int j = 0; j < num_samples; ++j) {
    for (int k = 0; k < n_vars; ++k) {
      if (l_lh_deterministic_test) {
        const int i_draw = (j + k + (iter - 1)) % 4;
        values(j, k) =
            (one / real_t(nt_repeat)) * (real_t(matrix(j, k)) + test_uniform_draws[i_draw]);
      } else
        values(j, k) = choose_permuted_random(nt_repeat, matrix(j, k), // In
                                              rng                      // Case state
        );
    }
  }
  Kokkos::deep_copy(X_u_one_lev, values);
}

//-----------------------------------------------------------------------
real_t choose_permuted_random(int nt_repeat, int p_matrix_element, // In
                              mt95& rng                            // Case state
) {
  // Draw one random jitter within the supplied permuted integer stratum.
  // nt_repeat: Number of samples before the sequence repeats
  // p_matrix_element: Permuted integer
  const real_t rand = rand_uniform_real(rng // Case state
  );
  return (one / real_t(nt_repeat)) * (real_t(p_matrix_element) + rand);
}

//-----------------------------------------------------------------------
void permute_height_time(int nt_repeat, int n_vars,         // In
                         SilhsInt2D one_height_time_matrix, // Out
                         mt95& rng                          // Case state
) {
  // Each variate receives an ordered draw of a new random permutation 0..N-1.
  // nt_repeat: Total number of sample points before sequence repeats.
  // n_vars: The number of variates in the uniform sample
  // one_height_time_matrix: nt_repeat x n_vars matrix of integers
  for (int i = 0; i < n_vars; ++i)
    rand_permute(nt_repeat, Kokkos::subview(one_height_time_matrix, Kokkos::ALL, i), // In / Out
                 rng                                                                 // Case state
    );
}

//-----------------------------------------------------------------------
void rand_permute(int n,            // In
                  SilhsInt1D pvect, // Out
                  mt95& rng         // Case state
) {
  // Art Owen, Quasi-Monte Carlo Sampling 1.3, following Devroye (1986).
  // No reseed: every swap consumes the next Native double-precision draw.
  // n: Number of elements to permute
  // pvect: Array of n numbers in random order
  auto host = Kokkos::create_mirror_view(pvect);
  for (int j = 1; j <= n; ++j)
    host(j - 1) = j;
  for (int j = n; j >= 2; --j) {
    const int temp = host(j - 1);
    double rand;
    rng.genrand_real3(rand);
    const int k = int(double(j) * rand + 1.0);
    host(j - 1) = host(k - 1);
    host(k - 1) = temp;
  }
  for (int j = 0; j < n; ++j)
    host(j) = host(j) - 1;
  Kokkos::deep_copy(pvect, host);
}
} // namespace clubb_kokkos
