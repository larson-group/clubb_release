// Implement CLUBB's named timing backends.
//
// timer_initialize, timer_disable, timer_start, timer_stop and timer_finalize
// mirror code_timer_module.F90. GPTL builds use its C API directly; other builds
// maintain a flat process-CPU-time registry with strict nesting and inclusive
// and exclusive totals. See code_timer_module.hpp for usage and output details.
#include "code_timer_module.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#ifdef GPTL
#include <gptl.h>
#else
#include "clubb_real.hpp"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <utility>
#include <vector>
#endif

namespace clubb_kokkos {
namespace {

// This explicitly permitted private state mirrors Native OpenMP threadprivate
// bookkeeping. Timer calls do not introduce Kokkos dependencies or GPU fences.
thread_local bool l_timer_initialized = false;
thread_local bool l_timer_enabled     = true;

#ifndef GPTL
constexpr int registry_chunk_size = 16; // Registry slots added at each growth. [-]
constexpr int stack_chunk_size    = 32; // Stack slots added at each growth. [-]
constexpr int display_name_length = 48; // Timer-name width in the CPU summary. [-]

struct timer_entry_t {
  std::string name;
  real_t      inclusive_time = real_t(0);
  real_t      exclusive_time = real_t(0);
  int         call_count     = 0;
};

struct timer_stack_entry_t {
  int timer_id      = -1;
  real_t start_time = real_t(0);
  real_t child_time = real_t(0);
};

// Vector sizes are Native allocated capacities, not logical entry counts.
// Registry IDs are zero-based; stack_depth and num_timers remain counts.
thread_local std::vector<timer_entry_t>       timer_entries;
thread_local std::vector<timer_stack_entry_t> timer_stack;
thread_local int                              num_timers  = 0;
thread_local int                              stack_depth = 0;

void        cpu_timer_start(const std::string& timer_name);             // In
void        cpu_timer_stop(const std::string& timer_name);              // In
int         get_timer_id(const std::string& timer_name, bool l_create); // In
void        ensure_registry_capacity(int required_capacity);            // In
void        ensure_stack_capacity(int required_capacity);               // In
void        print_cpu_timer_summary(std::ostream& timer_unit);          // InOut
std::string timer_display_name(const std::string& timer_name);          // In
void        require_empty_stack(const std::string& caller_name);        // In
void        reset_cpu_timer_state();
#endif

std::string       timer_absolute_path(const std::string& file_path); // In
[[noreturn]] void timer_fail(const std::string& message);            // In

} // namespace

// -----------------------------------------------------------------------------
void timer_initialize() {
  // Description:
  //   Initializes the active timing backend.

#ifdef GPTL
  int ret_code; // Return code from GPTL calls. [-]
#endif

  // --------------------------- Begin Code ---------------------------
  l_timer_enabled = true;

#ifdef GPTL
  if (l_timer_initialized) {
    ret_code = GPTLfinalize();
  }

  ret_code = GPTLsetoption(GPTLprint_method, GPTLfull_tree);
  ret_code = GPTLsetoption(GPTLabort_on_error, 1);
  ret_code = GPTLsetoption(GPTLoverhead, 0);
  ret_code = GPTLinitialize();
  if (ret_code != 0) {
    timer_fail("GPTLinitialize failed.");
  }
#else
  if (l_timer_initialized) {
    require_empty_stack("timer_initialize");
  }

  reset_cpu_timer_state();
#endif

  l_timer_initialized = true;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void timer_disable() {
  // Description:
  //   Disables the active timing backend and discards any recorded timer state.

#ifdef GPTL
  [[maybe_unused]] int ret_code; // Return code from GPTL calls. [-]
#endif

  // --------------------------- Begin Code ---------------------------
  if (l_timer_initialized) {
#ifdef GPTL
    ret_code = GPTLfinalize();
    (void)ret_code; // Native intentionally ignores this return code.
#else
    reset_cpu_timer_state();
#endif
  }

  l_timer_initialized = false;
  l_timer_enabled     = false;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void timer_start(const std::string& timer_name) { // In
  // Description:
  //   Starts or enters a named timer region. timer_name identifies the region. [-]

#ifdef GPTL
  [[maybe_unused]] int ret_code; // Return code from GPTL calls. [-]
#endif

  // --------------------------- Begin Code ---------------------------
  if (!l_timer_enabled) {
    return;
  }

  if (!l_timer_initialized) {
    timer_initialize();
  }

#ifdef GPTL
  // Fortran TRIM removes trailing blanks only; c_str supplies the C terminator.
  const std::string trimmed_name = timer_name.substr(0, timer_name.find_last_not_of(' ') + 1);
  ret_code                       = GPTLstart(trimmed_name.c_str());
  (void)ret_code; // Native intentionally ignores this return code.
#else
  cpu_timer_start(timer_name); // In
#endif
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void timer_stop(const std::string& timer_name) { // In
  // Description:
  //   Stops or exits the named region given by timer_name. [-]

#ifdef GPTL
  [[maybe_unused]] int ret_code; // Return code from GPTL calls. [-]
#endif

  // --------------------------- Begin Code ---------------------------
  if (!l_timer_enabled) {
    return;
  }

  if (!l_timer_initialized) {
    timer_fail("timer_stop called before timer_initialize or timer_start.");
  }

#ifdef GPTL
  const std::string trimmed_name = timer_name.substr(0, timer_name.find_last_not_of(' ') + 1);
  ret_code                       = GPTLstop(trimmed_name.c_str());
  (void)ret_code; // Native intentionally ignores this return code.
#else
  cpu_timer_stop(timer_name); // In
#endif
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void timer_finalize(const std::string& output_file) { // In
  // Description:
  //   Writes timer results to output_file and finalizes the active backend.

#ifdef GPTL
  [[maybe_unused]] int ret_code; // Return code from GPTL calls. [-]
#endif
  std::string output_path; // Absolute path used for timer output. [-]

  // --------------------------- Begin Code ---------------------------
  if (!l_timer_enabled) {
    return;
  }

  if (!l_timer_initialized) {
    return;
  }

  if (output_file.find_last_not_of(' ') == std::string::npos) {
    timer_fail("Timer output file must not be empty.");
  }

  output_path = timer_absolute_path(output_file); // In
  std::cout << "Writing timing results to: " << output_path << '\n';

#ifdef GPTL
  ret_code = GPTLpr_file(output_path.c_str());
  ret_code = GPTLfinalize();
  (void)ret_code; // Native intentionally ignores both return codes.
#else
  require_empty_stack("timer_finalize"); // In

  // A local stream replaces the Native NEWUNIT and STATUS="replace" file.
  std::ofstream file_unit(output_path, std::ios::out | std::ios::trunc);
  if (!file_unit) {
    timer_fail("Unable to open timer output file: " + output_path);
  }
  file_unit.exceptions(std::ios::badbit | std::ios::failbit);

  print_cpu_timer_summary(file_unit); // InOut
  file_unit.close();

  reset_cpu_timer_state();
#endif

  l_timer_initialized = false;
}
// -----------------------------------------------------------------------------

namespace {
#ifndef GPTL
// -----------------------------------------------------------------------------
void cpu_timer_start(const std::string& timer_name) { // In
  // Description:
  //   Starts or enters the cpu_time region identified by timer_name. [-]

  int    timer_id;   // Registry index for timer_name. [-]
  real_t start_time; // Process CPU time at timer start. [s]

  // --------------------------- Begin Code ---------------------------
  // std::clock measures process CPU time, as Native CPU_TIME does. Capture it
  // before registry/stack work. An unavailable process clock returns -1 second.
  const std::clock_t clock_value = std::clock();
  start_time = clock_value == std::clock_t(-1)
                   ? real_t(-1)
                   : static_cast<real_t>(clock_value) / static_cast<real_t>(CLOCKS_PER_SEC);

  timer_id = get_timer_id(timer_name, true); // In

  ensure_stack_capacity(stack_depth + 1); // In

  stack_depth = stack_depth + 1;
  timer_stack[stack_depth - 1].timer_id = timer_id;

  timer_entries[timer_id].call_count = timer_entries[timer_id].call_count + 1;

  timer_stack[stack_depth - 1].start_time = start_time;
  timer_stack[stack_depth - 1].child_time = real_t(0);
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void cpu_timer_stop(const std::string& timer_name) { // In
  // Description:
  //   Stops the named cpu_time region. It must match the most recent start.

  int    timer_id;  // Registry index for timer_name. [-]
  real_t stop_time; // Process CPU time at timer stop. [s]
  real_t elapsed;   // Elapsed inclusive timer duration. [s]

  // --------------------------- Begin Code ---------------------------
  const std::clock_t clock_value = std::clock();
  stop_time = clock_value == std::clock_t(-1)
                  ? real_t(-1)
                  : static_cast<real_t>(clock_value) / static_cast<real_t>(CLOCKS_PER_SEC);

  timer_id = get_timer_id(timer_name, false); // In

  if (stack_depth < 1) {
    timer_fail("timer_stop called with no active timers: " +
               timer_name.substr(0, timer_name.find_last_not_of(' ') + 1));
  }

  if (timer_stack[stack_depth - 1].timer_id != timer_id) {
    timer_fail("timer_stop nesting mismatch. Expected " +
               timer_entries[timer_stack[stack_depth - 1].timer_id].name + " but received " +
               timer_name.substr(0, timer_name.find_last_not_of(' ') + 1) + ".");
  }

  elapsed = stop_time - timer_stack[stack_depth - 1].start_time;

  timer_entries[timer_id].inclusive_time = timer_entries[timer_id].inclusive_time + elapsed;
  timer_entries[timer_id].exclusive_time =
      timer_entries[timer_id].exclusive_time + elapsed - timer_stack[stack_depth - 1].child_time;

  stack_depth = stack_depth - 1;

  if (stack_depth > 0) {
    timer_stack[stack_depth - 1].child_time = timer_stack[stack_depth - 1].child_time + elapsed;
  }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
int get_timer_id(const std::string& timer_name, bool l_create) { // In
  // Description:
  //   Returns the registry index for timer_name, optionally creating an entry.

  std::string trimmed_name; // Timer name without trailing blanks. [-]
  int timer_id;             // Registry index for timer_name, or -1 if not found. [-]

  // --------------------------- Begin Code ---------------------------
  trimmed_name = timer_name.substr(0, timer_name.find_last_not_of(' ') + 1);

  if (trimmed_name.empty()) {
    timer_fail("Timer name must not be empty.");
  }

  timer_id = -1;

  for (int i = 0; i < num_timers; ++i) {
    if (timer_entries[i].name == trimmed_name) {
      timer_id = i;
      return timer_id;
    }
  }

  if (!l_create) {
    timer_fail("Unknown timer name: " + trimmed_name);
  }

  ensure_registry_capacity(num_timers + 1); // In

  num_timers = num_timers + 1;
  timer_entries[num_timers - 1].name           = trimmed_name;
  timer_entries[num_timers - 1].inclusive_time = real_t(0);
  timer_entries[num_timers - 1].exclusive_time = real_t(0);
  timer_entries[num_timers - 1].call_count     = 0;

  timer_id = num_timers - 1;
  return timer_id;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void ensure_registry_capacity(int required_capacity) { // In
  // Description:
  //   Grows the timer registry to at least required_capacity slots. [-]

  std::vector<timer_entry_t> expanded_entries; // Temporary expanded registry. [-]
  int                        old_capacity;     // Existing registry capacity. [-]
  int                        new_capacity;     // Expanded registry capacity. [-]

  // --------------------------- Begin Code ---------------------------
  if (timer_entries.empty()) {
    timer_entries.resize(std::max(registry_chunk_size, required_capacity));
    return;
  }

  old_capacity = static_cast<int>(timer_entries.size());

  if (required_capacity <= old_capacity) {
    return;
  }

  new_capacity = std::max(required_capacity, old_capacity + registry_chunk_size);
  expanded_entries.resize(new_capacity);

  std::copy(timer_entries.begin(), timer_entries.end(), expanded_entries.begin());

  timer_entries = std::move(expanded_entries);
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void ensure_stack_capacity(int required_capacity) { // In
  // Description:
  //   Grows the active timer stack to at least required_capacity slots. [-]

  std::vector<timer_stack_entry_t> expanded_stack; // Temporary expanded stack. [-]
  int                              old_capacity;   // Existing stack capacity. [-]
  int                              new_capacity;   // Expanded stack capacity. [-]

  // --------------------------- Begin Code ---------------------------
  if (timer_stack.empty()) {
    timer_stack.resize(std::max(stack_chunk_size, required_capacity));
    return;
  }

  old_capacity = static_cast<int>(timer_stack.size());

  if (required_capacity <= old_capacity) {
    return;
  }

  new_capacity = std::max(required_capacity, old_capacity + stack_chunk_size);
  expanded_stack.resize(new_capacity);

  std::copy(timer_stack.begin(), timer_stack.end(), expanded_stack.begin());

  timer_stack = std::move(expanded_stack);
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void print_cpu_timer_summary(std::ostream& timer_unit) { // InOut
  // Description:
  //   Writes a flat summary of the cpu_time registry to the caller's stream.
  //   A local output stream replaces the Native integer timer_unit.

  real_t total_exclusive_time; // Sum of exclusive time over all timers. [s]
  real_t average_time;         // Inclusive time per timer call. [s]
  real_t exclusive_percent;    // Percentage of total exclusive time. [%]

  // --------------------------- Begin Code ---------------------------
  timer_unit.imbue(std::locale::classic());
  timer_unit << "========================= CLUBB TIMER SUMMARY =========================\n";
  timer_unit << "Backend: cpu_time\n";

  if (num_timers < 1) {
    timer_unit << "No timers were recorded.\n";
    timer_unit << "=======================================================================\n";
    return;
  }

  total_exclusive_time = real_t(0);
  for (int i = 0; i < num_timers; ++i) {
    total_exclusive_time += timer_entries[i].exclusive_time;
  }

  timer_unit << std::right << std::setw(48) << "Timer" << ' ' << std::setw(10) << "Calls" << ' '
             << std::setw(14) << "Inclusive(s)" << ' ' << std::setw(14) << "Exclusive(s)" << ' '
             << std::setw(14) << "Avg Inc(s)" << ' ' << std::setw(10) << "Excl Pct" << '\n';
  timer_unit << std::string(display_name_length, '-') << ' ' << std::string(10, '-') << ' '
             << std::string(14, '-') << ' ' << std::string(14, '-') << ' ' << std::string(14, '-')
             << ' ' << std::string(10, '-') << '\n';

  // C++ field widths are minima. Render each numeric field separately so an
  // over-width value becomes stars, as in Native I10/F14.6/F10.2 output.
  const auto field = [](const std::string& text, std::size_t width) {
    return text.size() > width ? std::string(width, '*')
                               : std::string(width - text.size(), ' ') + text;
  };
  const auto real_field = [&field](real_t value, std::size_t width, int precision) {
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::fixed << std::setprecision(precision) << value;
    return field(text.str(), width);
  };

  for (int i = 0; i < num_timers; ++i) {
    average_time =
        timer_entries[i].inclusive_time / static_cast<real_t>(timer_entries[i].call_count);

    if (total_exclusive_time > real_t(0)) {
      exclusive_percent = real_t(100) * timer_entries[i].exclusive_time / total_exclusive_time;
    } else {
      exclusive_percent = real_t(0);
    }

    timer_unit << timer_display_name(timer_entries[i].name) << ' '
               << field(std::to_string(timer_entries[i].call_count), 10) << ' '
               << real_field(timer_entries[i].inclusive_time, 14, 6) << ' '
               << real_field(timer_entries[i].exclusive_time, 14, 6) << ' '
               << real_field(average_time, 14, 6) << ' ' << real_field(exclusive_percent, 10, 2)
               << '\n';
  }

  timer_unit << "----------------------------------------------------------------------\n";
  timer_unit << "Total exclusive timed cpu_time seconds: "
             << real_field(total_exclusive_time, 14, 6) << '\n';
  timer_unit << "=======================================================================\n";
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
std::string timer_display_name(const std::string& timer_name) { // In
  // Description:
  //   Returns a fixed-width display name for the cpu_time summary. [-]

  std::string display_name(display_name_length, ' '); // Fixed-width display name. [-]

  // --------------------------- Begin Code ---------------------------
  if (timer_name.find_last_not_of(' ') + 1 <= display_name_length) {
    display_name.replace(0, timer_name.find_last_not_of(' ') + 1,
                         timer_name.substr(0, timer_name.find_last_not_of(' ') + 1));
  } else {
    display_name.replace(0, display_name_length - 3, timer_name.substr(0, display_name_length - 3));
    display_name.replace(display_name_length - 3, 3, "...");
  }

  return display_name;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void require_empty_stack(const std::string& caller_name) { // In
  // Description:
  //   Verifies that no cpu_time regions are active before caller_name proceeds.

  // --------------------------- Begin Code ---------------------------
  if (stack_depth > 0) {
    timer_fail(caller_name.substr(0, caller_name.find_last_not_of(' ') + 1) +
               " called while timer " + timer_entries[timer_stack[stack_depth - 1].timer_id].name +
               " is still active.");
  }
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
void reset_cpu_timer_state() {
  // Description:
  //   Clears the cpu_time registry and active stack, releasing their storage.

  // --------------------------- Begin Code ---------------------------
  timer_entries = std::vector<timer_entry_t>();
  timer_stack   = std::vector<timer_stack_entry_t>();

  num_timers  = 0;
  stack_depth = 0;
}
// -----------------------------------------------------------------------------
#endif

// -----------------------------------------------------------------------------
std::string timer_absolute_path(const std::string& file_path) { // In
  // Description:
  //   Returns the timer output path, using PWD for a relative file_path. [-]
  //   Preserve Native's 1024-character assignment before the caller's TRIM.

  std::string absolute_path; // Absolute timer output path. [-]
  std::string cwd;           // Current working directory from the environment. [-]
  int         env_status;    // Status equivalent to get_environment_variable. [-]

  // --------------------------- Begin Code ---------------------------
  if (file_path.front() == '/') {
    absolute_path = file_path.substr(0, file_path.find_last_not_of(' ') + 1);
  } else {
    const char* cwd_value = std::getenv("PWD");
    cwd                   = cwd_value ? cwd_value : "";
    env_status            = cwd_value ? (cwd.size() <= 1024 ? 0 : -1) : 1;

    if (env_status == 0 && cwd.find_last_not_of(' ') != std::string::npos) {
      absolute_path = cwd.substr(0, cwd.find_last_not_of(' ') + 1) + "/" +
                      file_path.substr(0, file_path.find_last_not_of(' ') + 1);
    } else {
      absolute_path = file_path.substr(0, file_path.find_last_not_of(' ') + 1);
    }
  }

  // C++ paths omit fixed-length blank padding; truncation still precedes TRIM.
  absolute_path = absolute_path.substr(0, 1024);
  absolute_path.erase(absolute_path.find_last_not_of(' ') + 1);
  return absolute_path;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
[[noreturn]] void timer_fail(const std::string& message) { // In
  // Description:
  //   Reports an unrecoverable timer error through the C++ entry's fatal handler.

  // --------------------------- Begin Code ---------------------------
  // The standalone std::exception handler prints this message and returns 1,
  // replacing Native stderr plus ERROR STOP without adding a second diagnostic.
  throw std::runtime_error("CLUBB timer error: " +
                           message.substr(0, message.find_last_not_of(' ') + 1));
}
// -----------------------------------------------------------------------------

} // namespace
} // namespace clubb_kokkos
