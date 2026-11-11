// Uniform Latin-hypercube samples using a caller-owned Native MT95 stream.
// Public functions: rand_uniform_real, generate_uniform_lh_sample,
// choose_permuted_random, permute_height_time and rand_permute.
#pragma once
#include "CLUBB_core/array_types.hpp"
#include "SILHS/mt95.hpp"
#include "SILHS/latin_hypercube_arrays.hpp"
namespace clubb_kokkos {
using SilhsInt1D = Kokkos::View<int*, Kokkos::LayoutLeft>;
using SilhsInt2D = Kokkos::View<int**, Kokkos::LayoutLeft>;
real_t rand_uniform_real(mt95& rng // Case state
);
void generate_uniform_lh_sample(int iter, int num_samples, int sequence_length, int n_vars, // In
                                bool l_lh_deterministic_test,                               // In
                                Real2D X_u_one_lev,                                         // Out
                                latin_hypercube_arrays_state& sampling_state, // Case state
                                mt95& rng                                     // Case state
);
real_t choose_permuted_random(int nt_repeat, int p_matrix_element, // In
                              mt95& rng                            // Case state
);
void permute_height_time(int nt_repeat, int n_vars,         // In
                         SilhsInt2D one_height_time_matrix, // Out
                         mt95& rng                          // Case state
);
void rand_permute(int n,            // In
                  SilhsInt1D pvect, // Out
                  mt95& rng         // Case state
);

template <class Field>
void generate_uniform_lh_sample(int iter, int num_samples, int sequence_length, int n_vars, // In
                                bool l_lh_deterministic_test,                               // In
                                Field X_u_one_lev,                                          // Out
                                latin_hypercube_arrays_state& sampling_state, // Case state
                                mt95& rng                                     // Case state
) {
  // Public strided/host output is staged by the shared array owner.
  // iter: Model iteration number
  // num_samples: `n'  Number of samples generated
  // sequence_length: `n_t' Num. random samples before sequence repeats
  // n_vars: Number of uniform variables to generate
  // l_lh_deterministic_test: Ordered strata and repeating draws for repeatable testing
  // X_u_one_lev: num_samples by n_vars matrix, X each row of which is a n_vars-dimensional sample
  auto execution = to_execution_space(X_u_one_lev, false);
  generate_uniform_lh_sample(iter, num_samples, sequence_length, n_vars, // In
                             l_lh_deterministic_test,                    // In
                             execution,                                  // Out
                             sampling_state, rng                         // Case state
  );
  from_execution_space(X_u_one_lev, execution);
}
template <class Field>
void permute_height_time(int nt_repeat, int n_vars, Field one_height_time_matrix, mt95& rng) {
  // nt_repeat: Total number of sample points before sequence repeats.
  // n_vars: The number of variates in the uniform sample
  // one_height_time_matrix: nt_repeat x n_vars matrix of integers
  auto execution = to_execution_space(one_height_time_matrix, false);
  permute_height_time(nt_repeat, n_vars, execution, // In / Out
                      rng                           // Case state
  );
  from_execution_space(one_height_time_matrix, execution);
}
template <class Field>
void rand_permute(int n,       // In
                  Field pvect, // Out
                  mt95& rng    // Case state
) {
  // n: Number of elements to permute
  // pvect: Array of n numbers in random order
  auto execution = to_execution_space(pvect, false);
  rand_permute(n, execution, // In / Out
               rng           // Case state
  );
  from_execution_space(pvect, execution);
}
} // namespace clubb_kokkos
