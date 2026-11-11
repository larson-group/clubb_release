# Independent ownership and driver cleanup review

Reviewed the live `src/driver_data.hpp`, `src/clubb_driver.hpp`, and
`src/clubb_driver.cpp` against the task's frozen driver sources under
`build/validation/driver_data_lifetime/baseline/`. No live edits or model builds
were made during this review. Native phase/callee access auditing belongs to
the separate driver source review; this review checks the storage adapter,
loop transformations, field mappings, unchanged code and documentation.

## Executable review

The adapter matches the independently tested draft header byte for byte.
Its explicit bind/sync/modify/handle/release contract is internally consistent.
It owns no global registry or model state. Dense mutable HostSpace fields alias
the existing allocation when the default execution space uses HostSpace;
strided fields use retained dense storage. For model real fields, the adapter's
`execution_view_type` is exactly the corresponding `Real1D` through `Real4D`
alias: both use `View<real_t... , LayoutLeft>` with the default execution space.
The original explicit-HostSpace View type need not itself be the same C++ type;
the converting View constructor shares its allocation. Actual Serial pointer
aliasing was checked in the focused ownership probe.

All 158 execution declarations, host/execution visitor entries and Native-named
advance handles match one for one in the same order, with no renamed mapping,
omitted handle or duplicate field. The local advance handles are ordinary Views;
the six new kernels do not capture the `ExecutionFields` owner or `this`.
The ice metadata indices are copied to scalar locals before the launch.

All six transformed loops enumerate columns first through
`i=index%extent(0), k=index/extent(0)`. Their iteration counts cover the same
logical cells as the previous column-inner loops. Thermodynamic and momentum
forcing retain expression association. Ice retains the conditional ice/snow/
graupel addition order and one-based metadata index conversion. Both cloud
number branches and the total-water statistics temporary retain their original
per-cell expressions. Two zeroing loops become full-array scalar deep copies;
they assign the same exact zero without changing surrounding branch/call order.

The proof script checks all six per-cell bodies at executable-token level after
only the stated index and metadata substitutions. After removing the reviewed
lifetime insertions and undoing the six kernel transformations, two zero fills,
temporary View type and unused local `k` deletion, the rest of the driver matches
the baseline's 35,384 executable tokens. Removing the new execution ownership
declarations and include leaves the header's original 2,943 executable tokens.
This proves unchanged calls, arguments, conditions and operation order outside
the reviewed adaptations; it does not prove whole-model numerical agreement.

Source snapshots, hashes, the normalization script and its passing result are
in `build/validation/driver_data_lifetime/owner_cleanup_review/`. The root agent
owns numerical validation and any later source changes.

## Documentation corrections sent to root

- The driver opening description still calls `get_stats` output a "native
  snapshot". It is now a C++ statistics snapshot.
- The driver opening claim that every View access converts retained Native
  loop indices needs qualification: the new execution kernels use zero-based
  indices directly.
- The core API opening description still says statistics retain f2k calls;
  that predates the completed C++ statistics port.

The new driver lifetime descriptions correctly explain host initialization,
persistent execution fields, reset transfers, bounded advances and cleanup.
No executable or capture correction is required by this review.
