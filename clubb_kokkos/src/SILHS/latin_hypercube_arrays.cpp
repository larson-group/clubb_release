// Release case-owned Latin-hypercube permutation storage.
#include "SILHS/latin_hypercube_arrays.hpp"
namespace clubb_kokkos {
void cleanup_latin_hypercube_arrays(latin_hypercube_arrays_state& sampling_state // Case state
) {
  // Native deallocates only this array; prior_iter is set by the next allocation.
  sampling_state.one_height_time_matrix = {};
}
} // namespace clubb_kokkos
