---
name: clubb-kokkos-maintenance
description: Maintain, build, compare or profile the existing CLUBB C++/Kokkos implementation and its runtime, numerical and support contracts.
---

# Maintain the CLUBB Kokkos implementation

Use when the checkout includes `clubb_kokkos/`, currently developed on the
`clubb_kokkos` branch. This guidance does not add that implementation to master.
For a checkout without it, use the generic C++/Kokkos language subskill below.

In a checkout with that backend, use this skill for repairs, build/runtime
questions, source synchronization, comparisons and approved performance work.
This skill owns implementation-specific decisions. Read the
`clubb_kokkos/README.md` for supported options, build/run
commands, current representation and checks. The
`clubb_kokkos/CORE_CONVERSION_PLAN.md` and
`clubb_kokkos/DRIVER_CONVERSION_PLAN.md` record scope,
architecture and historical acceptance; those results belong to their recorded
versions, not automatically to today's working tree.

For actual translation/expansion or a faithful source audit, also use the
[porting strategy](../clubb-porting/SKILL.md) and its fidelity checklist. Resume
the existing implementation instead of restarting completed milestones. Current
Fortran is its source and numerical authority; JAX is a secondary reference.
A different future port must name its own source/reference.

Generic indexing, layout, ownership, workspace, kernel mapping and performance
advice lives in the [C++/Kokkos language subskill](../clubb-porting/languages/cpp-kokkos/SKILL.md).
Read that material when relevant instead of duplicating it here. Maintenance
does not authorize unrelated cleanup, tuning, other backends, goals or external
writes. For documentation and test reporting use their existing repo owners.

## Current ownership and integration boundaries

- Mirror source names/tree and public argument correspondence. Preserve current
  zero-based C++ storage metadata and case/scoped workspace ownership. Source
  ordinals, serialized values and f2k wire fields convert at their established
  boundaries; do not reintroduce pervasive one-based computational storage.
- Production science, statistics and derived-type state are C++ owned. The
  default standalone does not require Native statistics/type services. Generic
  whole-procedure bridges have been retired; do not restore them as fallbacks.
- The optional `KOKKOS_BFB` matrix path retains the documented Fortran/LAPACK
  reference dependency. Inspect the README and CMake's requested/effective mode
  and warning behavior; do not assume the default custom factorization produces
  the reference bits or claim all configurations are computationally independent.
- Test-only fieldwise UDT converters remain under `clubb_kokkos/src/f2k/`.
  Reuse their generator/owner for relevant interface changes. Do not interpret
  their presence as permission for new production bridge calls. Follow the
  README's generation checks and preserve unrelated generated formatting.
- Use current shared grid, context, transfer and stats owners. Preserve deep
  snapshots versus aliases, active batch/mask/bounds semantics, stats cadence,
  diagnostics and repeated-case reset/cleanup. A local numerical-engine style
  waiver does not relax strict driver/public interfaces.
- Keep using the existing `run_scm.py -exe` route and shared case utilities.
  Broader upstream selector/build integration is separate scope, not a routine
  prerequisite for backend maintenance.

## Current comparison contract

Use the `clubb_kokkos/auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py`
and `clubb_kokkos/README.md`.
The regular gate is CPU Fortran versus CPU Kokkos, zero tolerances and strict
output comparison. Select the documented reference matrix mode when claiming
that BFB result. Do not inherit JAX's tolerances or curated windows.

Use a separate CPU execution-space build. A CUDA-enabled Kokkos installation
can retain CUDA as its default execution space even if Serial is also enabled;
the selected package/configuration, not an enabled-backend list, determines
which runtime executes. The README owns exact installed build examples.

Short multicolumn-statistics runs are the normal development cycle. Standard/all
relevant statistics and active nondefault/batch/scalar/forcing cases add coverage.
The README includes a 600-step performance-milestone example; this does not
override a task requiring full native duration. Record effective limits and
use the final duration/settings agreed for the actual work.

GPU comparisons are occasional milestones under this implementation's existing
policy: record scope and the newly integrated piece where saved GPU results
change, even slightly. GPU/CPU BFB is not the regular integration gate. Do not
trace GPU numerical causes unless requested. Device compilation/execution is
still required for device-relevant changes. A new port with another objective
must establish its own policy instead of inheriting this one.

The runner's defaults and locking are owned by the script; inspect them before
parallel launches. Distinct requested output roots do not remove its checkout
lock. Preserve strict field/schema/health checks and do not raise thresholds to
hide a known mismatch.

For failures use [numerical diagnosis](../clubb-numerical-diagnosis/SKILL.md).
Record deliberate forced association/FMA/library behavior in the
`clubb_kokkos/BFB_ADAPTATIONS.md`, with location,
ordinary expression, observed cause and before/after evidence. Keep pending
adaptations labeled pending. A proposed bit-changing optimization needs an
explained numerical cause and measured benefit before changing the contract.

## Performance and final evidence

Use optimized target-device builds and current profiling/timing utilities.
Bridge/test-only transfer timing is not production GPU-physics timing. Preserve
the distinction between requested payload counters, actual DMA traffic,
allocation counts and end-to-end throughput. Report workload, precision,
backend/device, initialization/stats/I/O scope, repetitions and contention.
Generic profiling/policy details stay with the language subskill.

Run relevant existing CTests/Python checks and integrated cases as documented;
add focused probes only for actual coverage gaps. For source-derived changes,
complete the [fidelity audit](../clubb-porting/fidelity/SKILL.md) separately from
numerical validation. Follow [testing/reporting](../clubb-testing-and-reporting/SKILL.md)
for exact version/case evidence and [documentation](../clubb-documentation/SKILL.md)
for comments and change messages. Do not call a historical checkpoint a test
of current source after later changes.
