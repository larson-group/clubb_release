// Check timer lifecycle and nesting edges outside the normal SCM call sequence.
// Both backends emit actual summaries; utilities/time_clubb.py is also exercised
// with real persistent-worker runs during integration validation.
#include "CLUBB_core/code_timer_module.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const std::string& message) {
  // Identify the lifecycle or output contract that failed.
  if (!condition)
    throw std::runtime_error(message);
}

template <class Operation> void expect_failure(Operation operation, const std::string& message) {
  // Only wrapper/CPU errors throw; GPTL's abort-on-error policy is unchanged.
  try {
    operation();
  } catch (const std::runtime_error& error) {
    require(std::string(error.what()).find("CLUBB timer error: " + message) == 0,
            "Timer failure lost its Native diagnostic");
    return;
  }
  throw std::runtime_error("Expected timer failure: " + message);
}

std::string read(const std::filesystem::path& file) {
  // Read the backend's real summary without substituting synthetic clock values.
  std::ifstream input(file);
  require(input.is_open(), "Timer output file was not created");
  return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}
} // namespace

int main(int argc, char** argv) {
  // Use one scratch directory for named output files. No Kokkos initialization
  // or GPU synchronization belongs in this host-only timer service.
  if (argc != 2)
    return 2;
  try {
    using namespace clubb_kokkos;
    const std::filesystem::path directory(argv[1]);
    std::filesystem::create_directories(directory / "path with spaces");
    const auto summary   = directory / "path with spaces" / "nested.timing";
    const auto fresh     = directory / "fresh.timing";
    const auto absent    = directory / "absent.timing";
    const auto discarded = directory / "discarded.timing";
    for (const auto& file : {summary, fresh, absent, discarded})
      std::filesystem::remove(file);

    timer_finalize(absent.string());
    require(!std::filesystem::exists(absent), "Uninitialized finalize wrote a summary");
    expect_failure([] { timer_stop("missing"); },
                   "timer_stop called before timer_initialize or timer_start.");
    timer_disable();
    timer_start("");
    timer_stop("");
    timer_finalize("");

    timer_initialize();
    timer_start("advance_clubb_to_end");
    for (int call = 0; call < 3; ++call) {
      timer_start("child   ");
      timer_stop("child");
    }
    // More than the Native 16 registry and 32 stack slots exercises both growth paths.
    for (int i = 0; i < 40; ++i)
      timer_start("grow_" + std::to_string(i));
    for (int i = 39; i >= 0; --i)
      timer_stop("grow_" + std::to_string(i));
    timer_stop("advance_clubb_to_end");
    timer_finalize(summary.string());
    const auto text = read(summary);
    require(std::regex_search(text, std::regex(R"((^|\n)[*!]?\s*child\s+3\s)")),
            "Nested calls or trailing-blank timer names changed");
    require(text.find("advance_clubb_to_end") != std::string::npos &&
                text.find("grow_39") != std::string::npos,
            "Registry or stack growth lost a timer");
#ifndef GPTL
    require(text.find("CLUBB TIMER SUMMARY") != std::string::npos &&
                text.find("Backend: cpu_time") != std::string::npos &&
                text.find("Inclusive(s)") != std::string::npos &&
                text.find("Exclusive(s)") != std::string::npos,
            "CPU summary no longer follows time_clubb.py's input format");
#else
    require(text.find("Stats for thread 0:") != std::string::npos &&
                text.find("Called") != std::string::npos &&
                text.find("Recurse") != std::string::npos && text.find("Wall") != std::string::npos,
            "GPTL output no longer follows time_clubb.py's input format");
#endif
    timer_finalize(absent.string());
    require(!std::filesystem::exists(absent), "Repeated finalize wrote stale timer state");

    timer_start("fresh"); // Start initializes a new backend after finalization.
    timer_stop("fresh");
    expect_failure([] { timer_finalize("   "); }, "Timer output file must not be empty.");
    timer_finalize(fresh.string());
    require(!std::regex_search(read(fresh), std::regex(R"((^|\n)[*!]?\s*child\s+\d+\s)")),
            "Timer state survived reinitialization");

#ifndef GPTL
    timer_initialize();
    expect_failure([&] { timer_finalize((fresh / "file.timing").string()); },
                   "Unable to open timer output file:");
    expect_failure([] { timer_start(""); }, "Timer name must not be empty.");
    expect_failure([] { timer_stop("unknown"); }, "Unknown timer name: unknown");
    timer_start("outer");
    timer_start("inner");
    expect_failure([] { timer_stop("outer"); }, "timer_stop nesting mismatch.");
    expect_failure([] { timer_initialize(); }, "timer_initialize called while timer inner");
    expect_failure([&] { timer_finalize(discarded.string()); },
                   "timer_finalize called while timer inner");
    timer_stop("inner");
    timer_stop("outer");
    expect_failure([] { timer_stop("outer"); }, "timer_stop called with no active timers:");
#else
    timer_initialize();
#endif
    timer_start("discard_me");
    timer_disable(); // Disable deliberately discards even an active region.
    timer_stop("discard_me");
    timer_finalize(discarded.string());
    require(!std::filesystem::exists(discarded), "Disabled timer wrote discarded state");
    std::cout << "Timer lifecycle, nesting, growth and output contracts pass\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
