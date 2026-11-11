// Own temporary storage for one C++ test after its Kokkos ScopeGuard.
//
// Construct ScratchWorkspaceOwner in main after Kokkos initialization and pass
// workspace explicitly to test helpers. The root scope retains returned Views.
// A single completion fence precedes root-scope/chunk destruction, including
// early returns and exception unwinding; no per-routine fence is added.
#pragma once
#include "scratch_workspace.hpp"
namespace clubb_kokkos::test {
struct ScratchWorkspaceOwner {
  ScratchWorkspace workspace;
  ScratchWorkspace::Scope lifetime;
  ScratchWorkspaceOwner() : lifetime(workspace.scope()) {
    // Reserve one test lifetime without allocating chunks on construction.
  }
  ~ScratchWorkspaceOwner() {
    // Complete queued work while the root scope, chunks and Kokkos remain alive.
    Kokkos::fence("CLUBB test scratch workspace completion");
  }
};
} // namespace clubb_kokkos::test
