---
name: clubb-porting
description: Plan, start, continue or synchronize a CLUBB port between languages using ordered milestones, source fidelity and integrated numerical comparisons.
---

# Port CLUBB in this order

Use this strategy for a new port, an expansion, source synchronization or a
faithfulness audit. Source and target can be any languages. Name the source
being translated and the trusted executable used as the numerical oracle;
they need not be the same implementation. Existing ports are useful secondary
references, not automatic authorities. Respect the requested discussion/edit
scope; this skill does not authorize a goal, budget or external actions.

The milestone order below is the default CLUBB migration strategy. It is not
an invitation to invent a fresh order. It resolves interoperability before
scaling and leaves statistics as a trusted comparison service until the end.
Change it only for a concrete project constraint discussed with the user.
Speed modes change scheduling within milestones, not their prerequisites.
For a core-only scope, stop after core acceptance and report the retained trusted
driver/statistics dependencies; do not expand into the remaining milestones or
claim full independence. Omitted milestones remain outside scope, not completed.

For an existing port, establish which milestones are already complete and
resume at the first relevant gap. Do not recreate a bridge, defer already
owned statistics or redo completed phases just to follow the list.

## Plan with the user before implementation

Read existing plans, current implementation maintenance guidance and recent
worklog entries. Resolve the following choices together, reusing decisions
already made in the conversation:

- Source revision, trusted reference/backend and target language/runtime.
- Scope: core only, required driver, remaining science/services, supported
  options, deliberately excluded packages and the meaning of full hookup.
- Objective: device execution, latency/throughput, gradients or another language
  capability, rapid completion, or a combination. Choose representative cases
  and the capabilities that must work during the port, not only at its end.
- Numerical contract: BFB where feasible, otherwise justified comparison
  criteria, duration and coverage. Read [validation](validation/SKILL.md)
  before establishing the baseline. Never infer tolerances from language names.
- Pause/review milestones and final validation. A milestone is a required
  checkpoint, not an automatic permission question; pause where agreed.
- Execution speed: ask for **serial**, **balanced (recommended)** or
  **fastest**, unless already chosen. Use balanced when the user leaves the
  choice to the agent or an unattended task needs a default; do not silently
  choose fastest. Read [execution](execution/SKILL.md) for the chosen mode.

Keep the agreed scope, milestone state, adaptations and evidence in the port's
plan and the repository worklog. Current support lists/build commands belong
to implementation maintenance, not language subskills.

## Milestone 1 — Runnable standalone calling the trusted backend

Build the smallest target-language executable/entry point that looks like
CLUBB's standalone to its caller and delegates real work to the trusted
backend. Check one case's inputs, output, failure status and cleanup. A shell
demonstration is a bootstrap, not the completed driver port.

Keep the new backend contained. Reuse `run_scripts/run_scm.py -exe` for a
compatible executable, shared namelist creation and existing output/comparison
utilities. For a different runtime, inspect the existing case-runner owner and
use a thin adapter. Do not copy the scripting stack or add backend conditionals
throughout callers. Broad build/selector integration waits until full hookup.

## Milestone 2 — All non-statistics types and real interoperability

Immediately settle every non-statistics UDT/object needed by the agreed driver
and physics scope, not just the first error object. Decide fields, precision,
index conventions, layout, allocation, ownership, optional values, inout/out
updates and case/process lifetime. Read the relevant source and target
[language subskills](#language-subskills) and [fidelity](fidelity/SKILL.md).

Prove real calls and round trips, including mutation, errors, array shapes and
reset/cleanup; use focused boundary checks where an ordinary model run cannot
expose a transfer defect. Without this foundation, every later driver call
creates another local converter, hidden state store or representation retrofit.

Put temporary bridges in a distinct owner with a strict contract. Convert at
the boundary, keep unported calls explicit and remove obsolete bridge paths as
target calls replace them. Avoid global scratch/storage imposed solely by a
binding tool; isolate and serialize any unavoidable binding-owned state.

**Do not port the statistics UDT/backend here.** Keep an opaque/reference-owned
handle and the minimal call/data interface needed to send statistics to the
trusted backend. Science should already produce comparable output without
reimplementing its observer. All other required types are settled now.

## Milestone 3 — Port the actual main driver

Translate the driver's initialization, advance/time-loop and cleanup bodies,
with calls to still-unported computations/services through the explicit bridge.
Keep source file/routine correspondence, call order, error handling and state
ownership. A thin whole-driver forwarder does not complete this milestone.

Establish the ordinary cycle: port, source-review, connect, compare and correct.
Prove a few meaningful replacements before scaling out in any speed mode.
For device/JIT/gradient objectives, test an actual target-owned piece early;
successful delegation to the reference backend proves no target capability.
Keep the normal execution path usable rather than building a second feature
demo driver that bypasses diagnostics or errors.

## Milestone 4 — Port the CLUBB core through the running call tree

Start from the core controller and replace reached routines in meaningful,
compatible groups. Keep lower unported calls explicitly delegated. Do not
port disconnected leaves or build a broad unit-oracle suite as prerequisites
for progressing through the model. Drafting other files can overlap once the
basic integration/comparison cycle is established.

Source-audit every promoted piece. Accept a compatible group when its integrated
comparison meets the agreed contract; a passing group needs no separate model
run per file. If it fails, narrow the additions and investigate the first cause.
Continue checking intended device execution and required language features as
the graph grows. Finish with a core source/cleanup audit and milestone comparison.

## Milestone 5 — Remaining driver science and non-statistics services

Port the agreed outer science packages and actual engine bodies, not merely
their wrappers. Then complete required input/readers, configuration, printing,
errors, timers and debug services. Input/IOP modules may contain science and
must not be dismissed as I/O. Preserve file/timing formats consumed by shared
tools. A local style waiver for an engine never spreads into strict interfaces.

Keep reference-backed statistics throughout. Remove each other obsolete whole
routine bridge as its real target owner is integrated. Retaining the trusted
statistics service keeps physics errors distinguishable from sampling/output
errors while these remaining computations change.

## Milestone 6 — Audit and stabilize before touching statistics

Confirm the agreed non-statistics science/services are complete, the supported
call graph reaches them, nondefault branches are accounted for, and the broader
integrated comparison passes. Inventory remaining reference calls and justify
each one. Review current source and tested versions separately. This checkpoint
prevents simultaneous physics and observer changes from obscuring failures.
Do not turn this checkpoint into an automatic commit or push.

## Milestone 7 — Statistics last, as one coherent replacement

Only after the preceding checkpoint, replace the statistics state, sampling,
accumulation, callers, lifecycle and output backend coherently. Preserve cadence,
budgets, masks, scalar/batch dimensions, reset/restart behavior, snapshot ownership
and the output schema as well as values. Plan compatible drafts if useful, but
avoid a long-lived mixture of competing statistics owners or partial formats.

Compare real output files against the trusted backend. Remove the statistics
bridge when the complete target backend passes. Statistics comes last because
it is the observation mechanism used to validate every earlier milestone;
changing that mechanism early makes a science mismatch harder to diagnose.

## Milestone 8 — Full hookup and agreed final validation

Close the source/branch inventory and remove obsolete migration machinery;
document any explicitly retained dependency. Verify the agreed case/options
matrix, effective durations, output schema, target-device functionality and
required language capabilities on identified current sources/executables.
Keep full-duration, short, edge, gradient and device evidence distinct.

If broader integration is in scope, connect through existing backend/runtime
owners now. Reuse namelist/run/diff logic and keep callers thin. Do not rebuild
the surrounding scripting system or expand host/CI work implicitly.

## Synchronization and performance follow-up

For source updates, inventory all relevant old/new deltas, including removals,
then trace signatures, types, callers, generated interfaces and diagnostics.
Default numerical agreement does not close this audit. Resume the existing
port's architecture; do not restart the initial migration sequence.

Choose sound target indices/layouts/loops/ownership from the beginning. Deeper
optimization needs the agreed objective and measured evidence. For an approved
set of easy improvements, visit each routine once and apply the whole agreed
set rather than repeatedly sweeping by category. This does not authorize
unrelated tuning. Speed settings never relax fidelity or numerical acceptance.

## Language subskills

Read only languages/runtimes relevant as source, target or bridge:

- [Fortran and C interoperability](languages/fortran/SKILL.md).
- [Python, NumPy and F2PY boundaries](languages/python/SKILL.md).
- [JAX transformations and functional state](languages/jax/SKILL.md).
- [C++ and Kokkos execution/storage](languages/cpp-kokkos/SKILL.md).

For an existing implementation, also read its maintenance owner:
[JAX](../clubb-jax-maintenance/SKILL.md) or
[Kokkos](../clubb-kokkos-maintenance/SKILL.md). These own current support,
build/runtime commands and implementation-specific checks. Generic language
advice does not override their actual public contracts.
