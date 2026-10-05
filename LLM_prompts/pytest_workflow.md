# Focused pytest coverage

Use this workflow when adding, changing, reviewing, pruning or wiring CLUBB
pytests into CI. Read [tests/README.md](../tests/README.md) for current suite
commands and Jenkins ownership. For real numerical comparisons, also use the
[numerical diagnosis workflow](numerical_regression_workflow.md).

## Choose the right test level

A pytest should check a specific behavior quickly: a numerical kernel against
an independent reference, an interface or serialization boundary, a failure
path, an invariant, or a small isolated runner/UI contract. Aim for seconds per
check and minutes per focused suite; measure expensive checks rather than
assuming that a short test file is cheap. JIT compilation can be appropriate
for a small kernel, but does not make a whole model run a unit test.

Anything that runs actual SCM cases or exercises a full workflow a user would
run from the command line belongs in higher-level `tests/`. Put component-specific
checks in that component's `tests/` folder; reserve repository-root `tests/` for
shared entry points, comparisons across implementations and general CLUBB
regressions. Give a workflow an explicit entry point, required inputs and
meaningful assertions. Preserve its assertions when moving it out of pytest.
Wiring a pytest suite does not authorize adding real-case or application/CLI
workflow stages to Jenkins; add those stages only within the requested CI scope,
otherwise document the manual validation route.

Starting a subprocess does not by itself require a higher-level test. A short,
bounded protocol regression, synthetic NetCDF file or browser-event simulation
can remain a pytest when it isolates a specific contract without running cases
or a full application workflow.

## Admission and ownership

- Put new agent-written checks under the owner's
  `pytests/auto_llm_generated_pytests/`. This includes new coverage added while
  editing an admitted module: put the new checks in a provisional module.
  Maintaining an existing admitted assertion does not require moving it.
- Use `tests/run_pytests.sh` as the pytest entry point. Its suite selection,
  environment preparation and generated-test admission are owned by the
  wrappers; do not add repository-root pytest configuration for these policies.
- Generated tests remain provisional. Wrapper defaults exclude them; use
  `tests/run_pytests.sh -SUITE -include_generated` to review them explicitly.
  A passing CI run does not promote a test. Human review promotes useful
  coverage by moving it into the parent `pytests/` or merging it into an
  existing admitted module, removing the provisional copy.
- Keep fixtures/helpers with their existing suite owner, outside generated
  test modules when shared. Use existing owner paths for source/fixture lookups
  so promotion into the parent folder preserves them. Do not import one test module from another or
  create a generic helper framework for one assertion.
- Every maintained pytest needs an explicit Jenkins owner, including a route
  for reviewing provisional tests. Use the ownership table in `tests/README.md`;
  update its existing pipeline and runner for a new suite. Do not leave useful
  checks discoverable only on an individual developer's machine.
- The initial API/JAX migration classifies all their existing tests as
  generated. Their Jenkins stages explicitly include them to preserve coverage
  while humans review them. This is an admission status, not proof of quality.

## Make failures informative

Assert behavior that could realistically break. Prefer independent numerical
oracles, boundary cases, error propagation, round trips through real interfaces,
and column independence. A finite derivative is a useful NaN guard, but does
not establish a correct derivative; use analytic values or finite differences
at smooth points when derivative correctness is the purpose. A copied formula
inside the test is not production coverage. Source-contract checks can protect
public fields/calls/porting requirements, but do not substitute for numerical
checks or the full faithful-port audit.

Avoid assertions that only restate constants, count incidental callbacks,
freeze cosmetic labels/layout numbers, or check that another test exists.
Testing a numerical validator with a representative fault is valuable when it
proves that the real validation detects that fault; repeated tests of those
meta-tests usually add little. Name the user-visible contract and why it matters.

Use `pytest.skip` for an expected optional dependency or capability, with a
specific reason. Never print "SKIP" and return, or swallow a loader/ABI error,
assertion, invalid call or missing checked-in source as a pass/skip. A required
CI build/input must fail setup if absent. Real-case oracles must generate or
require explicit inputs rather than quietly skipping an old local output path.

## Review and completion

Inventory pytest definitions and collected parameterized cases separately.
Trace what Jenkins actually invokes, including provisional inclusion, rather
than equating a filename or collection success with nightly execution. For a
removal, explain the obsolete contract or the surviving coverage; leave
standalone validation scripts and archived research source intact unless their
removal was requested. Compare failures with the starting revision in the
same environment before claiming a regression.

Run the affected suite in its owned environment and appropriate higher-level
workflow after a move. Verify both default exclusion and explicit inclusion of
provisional coverage. Report passes, skips, failures and unexecuted checks
separately, with timings and any known baseline failures. Validate changed
Jenkinsfiles and review discovery, path assumptions and the final diff.
