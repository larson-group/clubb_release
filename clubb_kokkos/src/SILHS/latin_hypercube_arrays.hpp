// Case-owned Latin-hypercube permutation and iteration diagnostic storage.
// Public function: cleanup_latin_hypercube_arrays releases the permutation.
#pragma once
#include "CLUBB_core/array_types.hpp"
namespace clubb_kokkos {
struct latin_hypercube_arrays_state {
  Kokkos::View<int**, Kokkos::LayoutLeft> one_height_time_matrix;
  int prior_iter{}; // Set on the first allocation, as in Native generator state.
};
void cleanup_latin_hypercube_arrays(latin_hypercube_arrays_state& sampling_state // Case state
);
} // namespace clubb_kokkos
