---
name: clubb-pytests
description: Add, audit, prune or wire CLUBB pytests into Jenkins, including focused checks, test placement and review status, and higher-level workflow ownership.
---

Use its linked owner documentation for current commands and CI coverage.

Respect the user's scope and existing authorization. Selecting this skill does
not turn a review into edits, approve external writes or create a goal.

# Focused pytest coverage

Use this workflow when adding, changing, reviewing, pruning or wiring CLUBB
pytests into CI. Read [human_reviewed_tests/README.md](../../../human_reviewed_tests/README.md) for current suite
commands and Jenkins ownership. For real numerical comparisons, also use the
[numerical diagnosis workflow](../clubb-numerical-diagnosis/SKILL.md).

## Choose the right test level

A pytest should check a specific behavior quickly: a numerical kernel against
an independent reference, an interface or serialization boundary, a failure
path, an invariant, or a small isolated runner/UI contract. Aim for seconds per
check and minutes per focused suite; measure expensive checks rather than
assuming that a short test file is cheap. JIT compilation can be appropriate
for a small kernel, but does not make a whole model run a unit test.

Anything that runs actual SCM cases or exercises a full workflow a user would
run from the command line is an integration test, not a pytest. Put
component-specific agent-written ones in that component's
`auto_llm_generated_tests/integration_tests/`; the repository-root
`human_reviewed_tests/` holds shared entry points, comparisons across
implementations and general CLUBB regressions. Give a workflow an explicit entry point, required inputs and
meaningful assertions. Its top-level description should explain the scenario
and assertions. For Fortran-derived checks, reference the source routines and
native test/harness being adapted, and identify any added backend-specific
assertions. Preserve its assertions when moving it out of pytest.
Wiring a pytest suite does not authorize adding real-case or application/CLI
workflow stages to Jenkins; add those stages only within the requested CI scope,
otherwise document the manual validation route.

Starting a subprocess does not by itself require a higher-level test. A short,
bounded protocol regression, synthetic NetCDF file or browser-event simulation
can remain a pytest when it isolates a specific contract without running cases
or a full application workflow.

## Placement, review status and ownership

- Tests are sorted by review status. Each component has
  `human_reviewed_tests/` (written or reviewed by a human) and
  `auto_llm_generated_tests/` (written by an agent, not yet reviewed), laid out
  as in [human_reviewed_tests/README.md](../../../human_reviewed_tests/README.md#test-layout).
  Everything under an `auto_llm_generated_tests/` folder is unreviewed,
  including the one inside the root `human_reviewed_tests/`.
- Put new agent-written tests under the owner's `auto_llm_generated_tests/`:
  pytests in `pytests/unit_tests/`, `clubb_jax` checks against Fortran through
  the Python API in `pytests/f2py_unit_tests/`, and real cases or whole-driver
  workflows in `integration_tests/`. This includes new coverage added while
  editing a human-reviewed test: put the new checks in a generated module.
- Only a human moves a test into `human_reviewed_tests/`, after reviewing it.
  Change human-reviewed tests only when the user asks. A passing CI run does
  not make a test reviewed.
- Use `human_reviewed_tests/run_pytests.sh` as the pytest entry point. Its suite
  selection and environment preparation are owned by the wrappers; do not add
  repository-root pytest configuration for these policies.
- Keep fixtures/helpers with their existing suite owner in
  `auto_llm_generated_tests/pytests/`, outside test modules when shared. Use
  existing owner paths for source/fixture lookups. Do not import one test
  module from another or create a generic helper framework for one assertion.
- Every maintained test needs an explicit Jenkins owner or a documented manual
  route. Each generated pytest suite has its own Jenkins stage, listed in the
  ownership table in `human_reviewed_tests/README.md`; identify the existing
  pipeline and runner that should own a new suite. Suggest any needed Jenkins
  integration; modify Jenkins test definitions only when the user requests that
  CI work, following [Jenkins test authoring](../clubb-jenkins/test-authoring/SKILL.md).
  Record pending integration so useful checks are not left discoverable only
  on one machine.

## Make failures informative

Assert behavior that could realistically break. Prefer independent numerical
oracles, boundary cases, error propagation, round trips through real interfaces,
and column independence. A finite derivative is a useful NaN guard, but does
not establish a correct derivative; use analytic values or finite differences
at smooth points when derivative correctness is the purpose. A copied formula
inside the test is not production coverage. Source-contract checks can protect
public fields/calls/porting requirements, but do not substitute for numerical
checks or the full faithful-port audit.

JAX pytests outside `f2py_unit_tests` folders must remain
independent of `clubb_python_api`, `clubb_python` and F2PY, including reference
oracles and runner preloads. Use analytic/NumPy references or supported JAX
contracts there; reading Fortran source text for a public schema is fine.
Routine-level comparisons against compiled Fortran belong in an
`f2py_unit_tests` folder (`-jax_f2py_unit_tests`; `-jax` ignores them), which calls
the `clubb_python` wrappers with the same inputs as the JAX routine and fails
setup without a Python API build. Whole-model comparisons use the
Fortran executable.
Numerical comparisons must assert the expected shape before comparing values
so accidental scalar broadcasting cannot pass as a full-array reference.

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

Use [testing and reporting](../clubb-testing-and-reporting/SKILL.md) to choose
relevant integrated checks and report their Jenkins links. Explain what focused
pytests checked and what confidence they add beyond the Jenkins result, rather
than reporting only a suite total or duplicating the same local run.

Inventory pytest definitions and collected parameterized cases separately.
Trace what Jenkins actually invokes, including which stages are commented
out, rather than equating a filename or collection success with nightly execution. For a
removal, explain the obsolete contract or the surviving coverage; leave
standalone validation scripts and archived research source intact unless their
removal was requested. Compare failures with the starting revision in the
same environment before claiming a regression.

Run the affected suite in its owned environment and appropriate higher-level
workflow after a move, and confirm it collects the same tests as before. Report passes, skips, failures and unexecuted checks
separately, with timings and any known baseline failures. Validate changed
Jenkinsfiles and review discovery, path assumptions and the final diff.
