---
name: clubb-porting-fidelity
description: Audit source-to-target routine, interface, comment and behavioral correspondence while allowing deliberate target-language representations.
---

# Source fidelity and target-language adaptations

Use with [the porting strategy](../SKILL.md) for implementation or an audit.
For a focused audit, apply this checklist directly without repeating completed
migration milestones. The source named in the port plan controls the mirror;
the numerical oracle controls the agreed comparison. A wrapper or a passing
default case does not complete a port.

## Before editing

Read the complete canonical source in scope and the target's relevant callers,
types and API owners. Record routines in order; signatures/intents/optionals;
major branches, calls and comments; computations; statistics/clipping/errors;
ownership and lifecycle. For substantial work keep this outline as an artifact.
Read active build/configuration gates before deciding a branch is inactive.
Existing ports are secondary references, not substitutes for reading the source.

Distinguish strict public/interface code from any explicitly relaxed engine.
A waiver is local to its approved scope and never waives numerical behavior.

## Preserve the correspondence

- Mirror source directories/files and routine names/inventory/order for the
  agreed surface, including private helpers. Avoid convenience splits, aliases
  and public-driver reordering. Respect target syntax constraints explicitly.
- Keep argument names, leading order, intent groups, required/optional status,
  outputs and inout ownership. Keep definition and caller contracts aligned;
  do not retain stale wrapper-only arguments or make required inputs optional.
- Prefer positional calls and corresponding multiline groups where safe;
  retain compact one-line calls. Necessary target syntax such as keyword-only
  parameters is an explained adaptation, not a reason to drift all callers.
- Preserve scientific dependency order, clipping/masking/conditional updates,
  initialization, diagnostic/statistics placement and fatal checks. Compare
  field-specific arguments individually; similar solve calls can differ.
  Remove duplicate work/statistics when ownership moved inside a callee.
- Preserve physical variable names, meaningful temporary values, kinds,
  formulas and rounding-sensitive associations. Avoid rename-only aliases,
  shape shortcuts and defensive extent changes that obscure source dimensions.
- Keep meaningful descriptions, references, warnings, comments/dividers and
  relative block locations. Use docstrings where comments cannot be copied
  literally. Preserve readable spacing and the source's meaningful grouping.
  Adapt literal index descriptions when the representation changes; retain
  their meaning rather than leaving a wrong comment in the target.
- Carry over beginning/end validation, debug-gated checks, detailed diagnostics,
  per-column versus whole-batch failure semantics, optional statistics and
  error returns. A reduced diagnostic/stub is incomplete, even if documented.
- Inspect nondefault and disabled branches, including scalar/species loops,
  perturbations, solver variants, clipping and statistics. Identify the actual
  gate for an excluded feature. Omit build-specific arguments only when outside
  the agreed surface; do not treat unknown macros as false.

## Choose target representations early

Faithfulness preserves the scientific/interface map; it does not require
carrying source-language inefficiencies into the target. Analyze what is cheap
or idiomatic in the source but expensive, unsafe or incompatible in the target.
Use the relevant language subskills before the type and driver milestones.

- Use target-native storage indices and constants, converting once at actual
  interop/serialization boundaries. Preserve counts, physical coordinates,
  option values, opaque IDs and meaningful sentinel semantics separately.
  Do not scatter source-index offsets throughout target computations.
- Choose array dimensions, strides/layout and iteration mapping together.
  A dimension permutation is allowed if logical field/axis correspondence and
  every transfer/caller remain explicit and correct; never reverse blindly.
- Express independent work with suitable array/vector/kernel primitives and
  true recurrences with ordered state. Preserve reduction/RNG dependencies and
  the numerical contract. Hoist/fuse work only when its dependency and option
  semantics are demonstrated, not merely because a default case still passes.
- Use target-appropriate immutable state or explicit owned context when needed;
  do not recreate unsafe globals just because the source used module state.
  Preserve source argument order and add only justified explicit lifecycle/state
  context at a consistent boundary, updating all callers and ownership records.
- Reuse case-owned persistent fields or scoped workspace when lifetimes permit.
  Preserve initialization, independent snapshots, aliasing and cleanup; a reused
  buffer must not expose stale data on rare/error paths.

Document an intended mapping once at its shared owner and explain non-obvious
adaptations near the boundary. These are supported design choices, not unfinished
TODOs. Use a local TODO only for actual incomplete work, with reason and removal
condition. Explain rounding-sensitive arithmetic choices beside the operation.

## Helpers, bridges and API owners

A target-only helper is justified by a real language construct (for example a
compiled loop body), ownership/interop boundary, deliberately reduced diagnostic
block, or a large genuinely identical repeated block. Keep its source mapping
visible. Do not hide field-specific calls, clipping, statistics or branch logic
behind generic convenience helpers. Remove one-line API forwarding, rename-only
unpacking, error/stat wrappers and repeated local dtype/layout converters unless
there is an actual boundary reason.

Use the established shared API/bridge, with one stable return contract. Repair
its missing arguments or outputs rather than accepting several possible result
shapes or bypassing it through generated bindings. Raw bindings belong only in
their designated wrapper owner. Normalize precision/layout at the boundary that
needs it; avoid bulk conversions on every computational routine's entry.

Temporary delegation is explicit in the migration inventory. Remove the old
whole-routine forwarder, imports and format adapters when the target body takes
ownership. Never add fallback/retry machinery to hide missing dependencies or a
broken target runtime. An already independent port does not regain a reference
runtime dependency merely because initial migrations use a bridge.

## Completion audit

After integration and numerical checks, re-read every changed target against
the full source outline. Check separately:

1. Routine/file inventory, order, definitions and every new helper.
2. Argument names/order/groups, optionality, outputs, inout and all callers.
3. Block/call order, field-specific arguments, comments/dividers and spacing.
4. Formula/kind/temporary, dimensions/index bounds, masks, initialization and
   recurrence/reduction/RNG semantics under the documented mappings.
5. Debug/error/diagnostic/statistics behavior and nondefault/build-gated branches.
6. Shared API contracts, raw-binding isolation, removed wrappers and obsolete
   fallbacks. Search for helpers/aliases/adapters and classify each occurrence;
   a search or parser is a prompt for review, not semantic proof.
7. Every exclusion, delegation or incomplete block has a concrete owner/gate and
   explanation; intended adaptations are distinguished from defects.
8. Actual runtime coverage and numerical/device/feature evidence are scoped to
   identified versions. Passing tests do not waive structural defects.

Report source/target revisions, reviewed routines/branches, justified adaptations,
commands/results, untested configurations and remaining work. Resolve in-scope
defects before claiming completion. Use [validation](../validation/SKILL.md)
and [testing/reporting](../../clubb-testing-and-reporting/SKILL.md) for evidence;
neither syntax checks nor a source audit alone establish numerical correctness.
