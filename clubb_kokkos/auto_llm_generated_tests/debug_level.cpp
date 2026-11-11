// Check signed debug thresholds against Native and isolate C++ host threads.
// Negative/extreme levels must retain the Native comparison, and each new host
// thread starts at zero without forwarding its setting to Native.
#include "CLUBB_core/error_code.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

extern "C" void native_debug_set_c(int level);
extern "C" int native_debug_check_c(int level);

namespace {
void require(bool condition, const char* message) {
  // Report the threshold or ownership contract that failed.
  if (!condition)
    throw std::runtime_error(message);
}
}

int main() {
  // No device execution is needed: both independent debug owners run on the host.
  try {
    using namespace clubb_kokkos;
    require(clubb_at_least_debug_level_api(0) && !clubb_at_least_debug_level_api(1),
            "C++ debug level does not start at zero");
    const std::array<int, 8> levels{-5, -1, 0, 1, 2, 3,
                                    std::numeric_limits<int>::min(),
                                    std::numeric_limits<int>::max()};
    for (int value : levels) {
      set_clubb_debug_level_api(value);
      native_debug_set_c(value);
      for (int requested : levels)
        require(clubb_at_least_debug_level_api(requested) == bool(native_debug_check_c(requested)),
                "C++ debug comparison differs from Native");
    }
    set_clubb_debug_level_api(2);
    bool thread_ok = false;
    std::thread other([&] {
      thread_ok = clubb_at_least_debug_level_api(0) && !clubb_at_least_debug_level_api(1);
      set_clubb_debug_level_api(7);
      thread_ok = thread_ok && clubb_at_least_debug_level_api(7) &&
                  !clubb_at_least_debug_level_api(8);
    });
    other.join();
    require(thread_ok && clubb_at_least_debug_level_api(2) &&
                !clubb_at_least_debug_level_api(3),
            "Debug metadata leaked between C++ host threads");

    native_debug_set_c(0);
    require(clubb_at_least_debug_level_api(2) && !native_debug_check_c(2),
            "C++ setter still forwards to the Native debug owner");
    set_clubb_debug_level_api(0);
    std::cout << "C++ debug thresholds and thread isolation pass\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
