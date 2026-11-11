# Persistent field transfer owner

The [HostDeviceField implementation](../../src/driver_data.hpp) keeps an existing
host field and a default-execution-space allocation together. The
[focused check](../../auto_llm_generated_tests/driver_data.cpp) is registered as
`clubb_kokkos_driver_data`. Original drafts are archived locally under
`build/validation/driver_data_lifetime/promoted_draft/`; the permanent files are
the implementation owners. Driver integration and model validation are separate
from the ownership probe.

This is a C++ ownership adaptation, with no numerical algorithm. Explicit
driver field lists remain in the driver. The helper does not choose model
phases, change arithmetic, track fields globally, or synchronize on handle
access. Existing host allocations remain usable by initialization and reset
routines. Dense CPU Views alias those allocations; other cases retain separate
execution storage. A host packing allocation is retained only when a device
transfer cannot directly use the host mapping. LayoutStride support covers
driver subviews; current `clubb_params` storage itself is dense.

## Interface and state

- `bind(host)` attaches the host allocation. An identical pointer, shape and
  strides preserve freshness. A changed host window becomes authoritative;
  execution storage is reused when extents match.
- `modify_host()` and `modify_device()` identify the latest writer. They do
  not copy or fence. Before partial or InOut writes, synchronize that side.
  A complete Out overwrite can deliberately discard the other side. There is
  no merge for conflicting independent writes.
- `sync_device()` and `sync_host()` transfer only if the other side is newer.
  They use blocking Kokkos deep copies and retain their allocations.
- `host_view()` and `execution_view()` return handles without synchronizing.
- `release()` clears both handles. Complete outstanding kernels before
  rebinding or releasing. Unbound access and synchronization throw
  `std::logic_error`.

New allocations retain normal Kokkos zero initialization. A persistent output
keeps unwritten entries on later calls; the driver audit must check conditionally
unwritten outputs against Native case lifetime. An identical host rebind does
not detect in-place resets: call `modify_host()` explicitly after reset.

## Focused checks

Both probes passed with NVHPC 26.9 and Kokkos 5.2.2: Serial and actual CUDA
execution on the RTX 3080 selected by its documented UUID. Checked rank 1–4 and
zero-sized fields, dense CPU aliasing, both transfer directions, strided window
updates without neighbor corruption, host modifications, same-size new windows,
stride-only rebinding, shape changes, and release. Kokkos Tools callbacks check
that already-synchronized calls perform no copy or allocation, and same-shape
window/stride changes allocate no storage. GPU values use exactly representable
integers, so these checks isolate transfers from math-library differences.

The isolated CMake project, binary/source hashes, commands and result scopes are
under `build/validation/driver_data_lifetime/owner_probe/RESULT.json`. No shared
build, live source or model harness was changed by this draft owner.
