// Centralized named timing calls for CLUBB.
//
// Usage:
//   timer_initialize() initializes the active backend before timed regions.
//   timer_disable() disables timing and discards recorded state.
//   timer_start(name) enters a region; timer_stop(name) exits that region.
//   timer_finalize(filename) writes timing results and finalizes the backend.
//
// Understanding:
//   GPTL builds forward to the GPTL C API. Other builds record process CPU time
//   in a dynamically growing flat registry. Regions must be explicitly nested:
//   if B starts inside A, B must stop before A. The CPU backend enforces this
//   order; names retain leading blanks and discard trailing blanks.
//
//   The CPU summary reports call count, inclusive time, exclusive time, average
//   inclusive time per call and exclusive percent. Inclusive time includes
//   nested regions; exclusive time subtracts them. GPTL controls its own output
//   format. Both outputs are consumable by utilities/time_clubb.py.
//
//   Private per-thread state mirrors Native OpenMP threadprivate bookkeeping.
//   Starting a timer initializes the backend lazily. Disabling discards active
//   regions; an explicit initialize reenables timing. Reinitialization and
//   finalization require an empty CPU stack. No device fences are added.
#pragma once

#include <string>

namespace clubb_kokkos {

void timer_initialize();
void timer_disable();

// timer_name: name of the timer region. [-]
void timer_start(const std::string& timer_name); // In
void timer_stop(const std::string& timer_name);  // In

// output_file: file used for timer output. Relative paths use PWD. [-]
void timer_finalize(const std::string& output_file); // In

} // namespace clubb_kokkos
