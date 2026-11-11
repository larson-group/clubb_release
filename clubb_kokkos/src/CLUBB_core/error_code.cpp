// Store the diagnostic threshold used by CLUBB's host-side checks.
//
// clubb_at_least_debug_level_api and set_clubb_debug_level_api mirror the Native
// accessors. Private per-thread diagnostic metadata follows Native threadprivate
// state, as with the timer service; no model fields are stored here.
#include "error_code.hpp"

namespace clubb_kokkos {
namespace {
thread_local int clubb_debug_level = 0; // Current host thread's diagnostic threshold [-].
}

bool clubb_at_least_debug_level_api(int level) { // In
  // Description:
  //   Checks to see if CLUBB has been set to the specified debug level.
  return level <= clubb_debug_level;
}

void set_clubb_debug_level_api(int level) { // In
  // Description:
  //   Set the diagnostic threshold. Native accepts any signed integer:
  //   0 prints no debug messages; 1 enables lightweight messages; 2 enables
  //   checks such as NaNs and spurious negative values. Do not clamp levels.
  clubb_debug_level = level;
}

} // namespace clubb_kokkos
