---
name: clubb-porting-execution
description: Select and run serial, balanced or fastest porting pipelines with compatible batches, bounded workers and isolated validation.
---

# Execution speed

Use with [the ordered milestones](../SKILL.md). Ask once during planning unless
the user has already chosen; recommend balanced and use it when a default is
needed. Do not infer fastest from urgency alone. All modes establish the small
interoperability/driver/comparison flow before scaling, preserve source review
and use the same numerical and final acceptance requirements.

## Serial

One agent ports a meaningful piece, reviews its source correspondence, connects
it, tests the integrated model, corrects failures and repeats. Do not draft the
next piece while that gate is pending. Compiler-internal parallelism and the
target model's scientific parallelism are separate from this scheduling choice.
Use small development gates and the agreed larger milestone checks.

## Balanced — recommended default

One agent starts a frozen integration/test gate, then drafts or audits the next
independent pieces outside its source tree. Integrate compatible reviewed pieces
together and compare the group. Accept a passing group without per-file model
gates; narrow a failing group before diagnosing its first difference.

Overlap useful preparation and independent comparisons when supported safely
by the maintained runner. Revisit pending gates before relying on their results.
No subagents are implied by this mode. Avoid edits to files an active build uses;
use a frozen lane or keep new drafts separate.

## Fastest

Use a coordinator and a bounded pool of subagents once the initial pipeline
works. Assign disjoint complete source owners/routine families with source
revision, allowed files, shared contracts and required fidelity. The coordinator
owns live integration, shared build/generator/API changes and the comparison
schedule. Coordinate types/global indices/caller closures before dividing work.

Draft workers should read the full source, implement the agreed surface and do
minimal useful parse/import/compile checks where inexpensive. Add a tiny probe
only for an edge/option/ownership risk integrated cases would miss, or when
needed to be reasonably confident integration will expose a defect. Do not
build a private numerical harness or prove every routine in isolation before
handoff. Source review and honest limitations still matter.

The coordinator should:

- Keep reviewed drafts ready while comparisons run, integrate compatible groups
  and run multiple independent test tasks/lanes when resources and runner locks
  allow. Use the existing harness, not a new per-worker runner.
- On failure, narrow the additions, find the first divergence and fix small
  clear defects promptly. Send substantial repairs back to the owner with a
  reproducer, exact versions/settings, failing operands/contract and allowed
  files; continue independent useful work while it is repaired.
- Check returned deltas against current live callers and recheck combined code.
  A worker's standalone pass is not combined-model evidence.
- Reuse idle workers with explicit reassignment. Nested delegation consumes the
  same slot/resource budget; it is not unlimited additional capacity.

If delegation is unavailable or not permitted, report that constraint and use
balanced scheduling; never claim subagent review occurred. Do not assign fixed
model names or worker counts in reusable guidance.

## Contracts that make overlap reliable

Keep agent topology, build/test process workers, model batch width and target
kernel parallelism distinct. Bound concurrency by CPU, RAM/compiler memory, GPU
memory/contention and disk; more simultaneous jobs can reduce throughput.

For each handoff record owned files/base and hashes, routine inventory,
adaptations/dependencies, checks actually run, review provenance, limitations and
integration steps still owned by the coordinator. Distinguish drafted, audited,
compiled, integrated, exercised and numerically accepted. An author's second
pass is not an independent review. If independent review is required, use an
actually distinct reviewer; an old worker name does not establish independence.

Use isolated source, build, runtime-code and output owners for simultaneous
validation generations. Inspect runner locks: separate output names alone may
not permit concurrency if the harness locks the checkout. Do not bypass those
locks. Share immutable input datasets with preserved symlinks when appropriate;
do not accidentally dereference a large benchmark tree into every lane.

Freeze and identify each tested source/binary/input/environment pair. Review
minimal passing deltas into the live tree instead of copying a lane over it,
then recheck the combined model. Later edits require relevant fresh evidence or
an explicit scoped compatibility argument; never relabel old tests as current.

Plan compact manifests/log retention and disk use. Keep failed setup/config
attempts separate from science failures and respect authorization for cleanup.
Finalize on a source freeze with the required completed checks; a queued or
running comparison is pending evidence.
