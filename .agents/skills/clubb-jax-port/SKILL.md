---
name: clubb-jax-port
description: Port or audit Fortran-derived JAX code in CLUBB, including driver integration and structural fidelity.
---

Its referenced standards are part of the workflow when applicable.

Respect the user's current scope and prior authorization. A discussion or
review request remains a discussion or review; selecting this skill does not
approve broader changes, external writes, or creating a goal.

# Port or audit JAX code

Use for translating Fortran into `clubb_jax/src`, wiring a port into the JAX
driver, or auditing such a port. Use
[the faithful-port standards](../clubb-fortran-port/SKILL.md)
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
[clubb_jax/README.md](../../../clubb_jax/README.md) for normal runtime use and
[clubb_jax/JAX_CONVERSION_PLAN.md](../../../clubb_jax/JAX_CONVERSION_PLAN.md) for
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
  `clubb_jax/human_reviewed_tests/test_file_name_parity.py` enforces the
  file-for-file mirror; it is human-maintained ground truth. When you add,
  move, rename or remove a JAX source file, run
  `.venv-jax/bin/python -m pytest clubb_jax/human_reviewed_tests/test_file_name_parity.py`
  (after `python3 clubb_jax/run_jax.py -init_env`) before and after your
  change. Fix any new failure caused
  by a rename, merge or split you introduced by changing the JAX files. Do not
  edit the test, including its allowlists. If your change makes an allowlist
  entry stale or needs a new exception, stop and report the test output to the
  user so a human can update the test or revise the change.
- Generate the docstrings of ported routines with
  [port_fortran_declarations_to_jax_docstrings.py](../../../utilities/port_fortran_declarations_to_jax_docstrings/README.md)
  after adding or changing a JAX routine with a Fortran counterpart, or its
  Fortran comments or arguments:
  `python3 utilities/port_fortran_declarations_to_jax_docstrings/port_fortran_declarations_to_jax_docstrings.py -write_description -modules <module>`.
  It writes the full Fortran description, references, returned values and
  argument declarations. Do not edit that text in the JAX file; record
  adaptations of Fortran text that does not apply to JAX (Fortran-only names
  such as LAPACK, unported branches, other real JAX differences) in
  `utilities/port_fortran_declarations_to_jax_docstrings/custom_ports_of_fortran_descriptions.jsonc`
  and rerun. Fix wrong or nonstandard Fortran comments in the Fortran instead.
  Its `-check` must pass; the clubb_jax Jenkins job runs it.
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

## CLUBB coding standards in Python

These rules are the Python form of the
[CLUBB coding standards](../../../doc/clubb_coding_standards.md) that matter in
JAX code; numbers refer to that document. Its other standards are not
requirements for Python code. These rules do not apply to external code such as
`clubb_jax/src/Microphys/Morrison_microphys/` and
`clubb_jax/src/Radiation/BUGSrad/`.

Apply them to the code you add or change, not as a reason to reformat or clean
up other code. In ported code, the faithful-port rules win: do not rename,
reorder or restructure against the Fortran to satisfy a rule, and keep a Fortran
`known magic number` or `known magic flag` comment with its value. Follow these
rules where the Fortran does not decide, such as in JAX-only code and new
adaptations.

Names:

- Start logical flags with `l_` (#12, #34).
- Give a variable its base name on its home levels. Add `_zm` or `_zt` when it
  is moved to the other levels and on variables with no home level (#37).
- Use `k` for vertical, `i` for grid column, `j` for scalar, `sample` for
  subcolumn and `p` for PDF dimension indices (#44).

Values:

- Do not add unnamed numbers (#18). Import shared values from
  `clubb_jax/src/CLUBB_core/constants_clubb.py` instead of redefining them, so
  that a value cannot differ between files. Mark an unavoidable value, such as
  one from a case specification, with `# known magic number`.
- Do not pass a bare `True` or `False` to a CLUBB routine (#19). Assign a named
  `l_` variable just above the call, as the Fortran usually does, or mark an
  unavoidable value with `# known magic flag`.
- For dtypes and fixed precision (#39), follow the description in
  `clubb_jax/src/CLUBB_core/clubb_precision.py`.

Data flow:

- Avoid new module variables (#5, #24). A port may mirror a Fortran module
  variable, but traced code must not change one. A value read during tracing
  is frozen into the compiled function, so read only static configuration,
  such as the debug level, that way.
- Pass individual flags rather than `clubb_config_flags` or
  `silhs_config_flags`, unless the Fortran routine takes the whole type (#42).
- Keep each stats or file write immediately after its value is set or at the
  end of `advance_clubb_core`, so that restarts stay bit-for-bit (#43).
- Assign each output to a named variable and return that variable, rather than
  returning a calculation such as `return x * x + 5`. In a port, use the
  Fortran output name: the output argument or the function result name.
  Returning existing variables, a tuple of them or a plain call is fine, and
  short functions nested inside another function, such as `lax.cond` branches
  and loop bodies, are exempt.

Expressions and output:

- `jnp.where` evaluates both branches, unlike a Fortran `if` or `where` (#41).
  An unselected branch that divides by zero or takes a power of zero or a
  negative number can still produce NaN or infinite gradients. Make its input
  safe with an inner `jnp.where`, as `_pos_pow` does in
  `clubb_jax/src/Microphys/KK_microphys/PDF_integrals_means.py`.
- Put non-fatal messages inside a `clubb_at_least_debug_level` check, and write
  fatal messages to `sys.stderr` (#35). A plain `print` in traced code runs only
  during tracing, so print traced values with `jax.debug.callback`.

## Verification and completion

Judge JAX numerical correctness against the corresponding Fortran implementation;
preserving previous JAX output is not currently its acceptance standard. Use
[testing and reporting](../clubb-testing-and-reporting/SKILL.md#choose-checks-that-exercise-the-change)
for comparison evidence and the canonical [bit-change policy](../clubb-bit-changing/SKILL.md)
for declarations.

For comparison failures, use [numerical diagnosis](../clubb-numerical-diagnosis/SKILL.md)
to establish the reference and first divergence before changing acceptance.

Use the managed JAX launcher/environment documented in the README, rather than
hard-coded paths to another checkout's venv. Begin with imports/focused checks
and short comparisons while debugging. Keep cases/settings equal in JAX and
Fortran; record overrides. Never relax tolerances, disable stats needed for a
comparison, shorten the final requested run, or change physics to obtain a pass
without making the limitation explicit and getting acceptance where required.

For example, from the repo root:

```bash
./human_reviewed_tests/run_jax_vs_fortran_cases.py -cases bomex -workers 1 -max_iters 20
./human_reviewed_tests/run_jax_vs_fortran_cases.py -cases bomex -workers 1
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

Use [testing and reporting](../clubb-testing-and-reporting/SKILL.md) for final
check selection and reporting Jenkins links, extra local coverage and failures.

Finish with the faithful-port completion audit, including comments, spacing,
argument grouping, every new helper, no new failures of the filename-parity test
and a passing `port_fortran_declarations_to_jax_docstrings.py -check`.
Check the code you added or changed against the "CLUBB coding standards in
Python" section above. Report structural and numerical results separately, remaining support gates
and untested configurations. Do not declare
complete while known in-scope structural defects remain. Do not redesign output,
optimize hardware performance, or refactor unrelated files just to polish a port.

## Scope boundaries

- For JAX work, the JAX workflow is primary and the faithful-port standards
  supply source matching. Existing/external JAX files are secondary references;
  selecting this skill does not create a file-level goal or require cloning.
