# Port or audit JAX code

Use for translating Fortran into `clubb_jax/src`, wiring a port into the JAX
driver, or auditing such a port. Use
[the faithful-port standards](port_underlying_fortran_to_other_languages.md)
for routine names/order, signatures, calls, comments, formatting, branches and
acceptance. Existing JAX or external CLUBB-JAX files are secondary references;
current source Fortran is the authority.

## Scope and source outline

Infer the requested files, integration path and validation scope from the
current conversation. Ask only about unresolved choices that materially change
the work. A discussion-only request remains discussion. Do not create a goal,
set a token budget, clone an external repo, or expand into all mirrors merely
because this workflow is selected.

Before editing, outline the full Fortran routines, signatures, call ordering,
branches, stats/error behavior and ownership of returned/inout data. Include
module lifecycle and callers when changing the driver. Read
[clubb_jax/README.md](../clubb_jax/README.md) for normal runtime use and
[clubb_jax/JAX_CONVERSION_PLAN.md](../clubb_jax/JAX_CONVERSION_PLAN.md) for
whole-core architecture and acceptance. Confirm support claims against current
initialization gates and tests; file-specific historical/ignored notes are
optional examples, not a competing current plan.

Prefer one integrated call-tree branch at a time: replace a meaningful node,
validate the runnable model path, resolve that branch's shape/tracing/compilation
issues for the intended support level, then reuse its leaves elsewhere. Avoid
breadth-first disconnected ports. In mature JAX-owned standalone code, this
method does not recreate the old transitional Fortran delegation.

Identify strict interface/public-contract files and any explicitly relaxed
numerical cores at the start. Retain the public file/routine/argument/call
surface required for those cores; local helper freedom does not spread into
their strict interface or waive numerical coverage.

## JAX implementation boundaries

- Mirror source-tree/file layout under `clubb_jax/src/`; keep routine names,
  order, arguments, meaningful comments, dividers and readable spacing.
- Use current JAX state/types/grid infrastructure and same-named source types.
  Reuse shared infrastructure rather than creating per-file conversion layers.
- Keep physics array math in JAX. Use fixed shapes, functional updates, masks
  and appropriate JAX control flow. Preserve source semantics and BFB-sensitive
  ordering; explain necessary indexing, batching or immutable-state adaptations.
- Inspect non-default, debug, disabled and error paths rather than assuming the
  smoke case covers them. Preserve supported behavior and document unsupported
  branches next to their explicit gates; avoid silent stubs or fallback paths.
- The supported standalone is JAX-owned and requires no compiled Fortran library.
  Do not introduce a Fortran/F2PY fallback into that path. The original Fortran
  executable is the numerical comparison reference. A task explicitly scoped to
  transitional code must name any remaining API call and why it remains; that
  historical allowance does not authorize a new dependency in standalone code.
- Do not import raw generated bindings into physics ports. If Python API work
  is explicitly in scope, repair its public contract instead of hiding its
  dtype/layout or multiple-return-shape problems in a local adapter.
- Remove swappable driver fallbacks for a layer that is now JAX-owned. Confirm
  the active standalone/core call path actually reaches the new port.

## Verification and completion

For comparison failures, use [numerical diagnosis](numerical_regression_workflow.md)
to establish the reference and first divergence before changing acceptance.

Use the managed JAX launcher/environment documented in the README, rather than
hard-coded paths to another checkout's venv. Begin with imports/focused checks
and short comparisons while debugging. Keep cases/settings equal in JAX and
Fortran; record overrides. Never relax tolerances, disable stats needed for a
comparison, shorten the final requested run, or change physics to obtain a pass
without making the limitation explicit and getting acceptance where required.

For example, from the repo root:

```bash
./tests/run_jax_vs_fortran_cases.py -cases bomex -workers 1 -max_iters 20
./tests/run_jax_vs_fortran_cases.py -cases bomex -workers 1
```

The second command removes only the iterative cap; the curated case itself may
have configured limits. Inspect the harness's effective settings before calling
it a native-duration run. Use its default bindiff threshold unless the user
requests another. GPU comparisons use one worker. Final verification should
match the scope agreed for the file/driver, with native-case or broader checks
when that behavior is part of the claim. Focused tests alone do not establish
full driver parity, JIT support, or differentiability. For parameter/tuner work,
check runtime-value changes reuse compilation rather than making each parameter
set a new static specialization. Shape/dtype/physics changes may legitimately
specialize; measure reuse rather than infer it from a decorator.

For differentiation, use the normal driver and current README-supported static
debug/stats configuration instead of a second execution path that bypasses
diagnostics. Preserve established negative-debug continuation at the relevant
result/loss boundary. Preserve per-column failure status and healthy batch results
in tuner work; initialization/configuration errors still need their existing
stops. Repair failure semantics before adding retries. Keep ordinary tests and
comparisons loud. A short gradient
probe is not evidence of long-run tuning or all-physics differentiability.
For approved stochastic adaptations, distinguish algorithm/permutation parity
from random-stream identity and record the deterministic test mode and limits.

Finish with the faithful-port completion audit, including comments, spacing,
argument grouping and every new helper. Report structural and numerical results
separately, remaining support gates and untested configurations. Do not declare
complete while known in-scope structural defects remain. Do not redesign output,
optimize hardware performance, or refactor unrelated files just to polish a port.
