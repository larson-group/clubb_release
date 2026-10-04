# Diagnose numerical comparisons

Use for JAX/Fortran, compiler, BFB or multi-column comparison failures. A smoke
pass and a longer failure describe different tested behavior; neither establishes
whether the current change caused the failure.

## Establish a comparable reference

Record source revisions, build freshness, compiler/optimization flags, dtype,
backend, cases, timestep/duration, overrides, random settings, stats and bindiff
criteria. Rebuild the Fortran reference after source changes; do not compare a
new port with a stale executable. Run the starting revision with equivalent
settings when attribution matters. Check the harness's effective limits rather
than assuming an uncapped command means a native-duration case.

Prefer live comparisons with the current Fortran source as the numerical
oracle. Keep focused edge/contract tests for behavior that integrated cases do
not exercise, rather than freezing another large output fixture or testing a
copy of the implementation. State active steps/fields/species and major omissions;
a quiet or reduced-stats run provides limited coverage.

When changing a comparison harness, use fresh isolated outputs, verify the
unmodified control passes, and exercise an active representative defect that
causes the intended numerical/schema failure while both model runs complete.
A crash, stale file, metadata-only mismatch or unexercised mutation does not
prove numerical detection. Keep the documented default schema warnings, strict-mode failures and empty-
comparison safeguards aligned with the owner contract. For a regression guard,
confirm it detects the buggy baseline; a test-only improvement may be enough.

## Locate the first cause

Use short runs to find the first diverging timestep, field and routine, then
trace operands and call/return/state ownership to the source. Inspect disabled,
error, diagnostic and non-default branches as well as the default case. Use
cases with nonzero activity in the feature under investigation: a quiet cloud
or radiation path can leave an incorrect implementation unexercised.

For suspected roundoff, compare controlled optimization/precision configurations
and examine absolute/relative/ULP behavior before attributing accumulated drift
to harmless noise. Preserve ordering during diagnosis. Ensembles, perturbations
or relaxed thresholds can be useful agreed experiments; selected per-case
tolerance/duration exceptions remain local to their case definitions. Check that a diagnostic
would still detect a deliberately wrong result before trusting it.

## Preserve the meaning of a pass

Do not silently shorten the final run, loosen tolerances, disable relevant
stats or change physics to turn a failure green. An approved diagnostic is not
a permanent acceptance change. Keep exploratory overrides separate from the
requested validation, and record both. Report numerical and structural findings
separately, with the tested scope and remaining uncertainty.

For performance investigations, retain a compact reproducible timing/profile
summary plus separate detailed logs, including configuration and repetitions.
Reuse the existing timing/profile utilities; avoid a new private instrumentation
path in every kernel. Performance evidence does not waive source fidelity.
