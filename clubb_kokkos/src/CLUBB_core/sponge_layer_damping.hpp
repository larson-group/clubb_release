// Retain the Native sponge namelist value type for standalone configuration.
//
// The driver owns six instances, initializes their Native defaults and rejects
// enabled sponge damping, as in the supported JAX/Kokkos option set. Damping
// profiles and numerical routines remain outside this port's supported scope.
// The C++ namelist reader preserves the original derived-member input syntax.
#pragma once
#include "clubb_real.hpp"

namespace clubb_kokkos {

struct sponge_damp_settings {
  real_t tau_sponge_damp_min{}; // Minimum damping time at model top [s].
  real_t tau_sponge_damp_max{}; // Maximum damping time at layer base [s].
  real_t sponge_damp_depth{};   // Depth as a fraction of domain height [-].
  bool l_sponge_damping = false; // Enable sponge damping [-].
};

} // namespace clubb_kokkos
