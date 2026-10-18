---
name: clubb-bit-changing
description: Decide whether a Fortran answer change needs a BIT_CHANGING declaration, identify affected BFB tests or option conditions, and separate it from result-preserving commits.
---

# Declare bit-changing commits sparingly

Use before splitting commits or writing/revising a bit-change declaration.
This is the canonical policy for `BIT_CHANGING`; message-writing, testing and
finalization workflows link here.

## Purpose and scope

The marker records intentional numerical answer changes so that future changes
in Fortran host models, especially CAM, can be traced to CLUBB commits. These
changes usually originate in `src/CLUBB_core`. Decide from the calculation and
its consumers, rather than the directory name alone.

Fortran standalone answer changes can also need the marker: changing prescribed
forcings may break CLUBB's BFB tests while leaving host answers unchanged. State
that standalone scope explicitly. The [discussion on PR #1399](https://github.com/larson-group/clubb/pull/1399#discussion_r4303880250)
retains this exception.

JAX-only answer changes currently do not receive the marker. JAX correctness
targets agreement with the corresponding Fortran implementation; preserving
earlier JAX answers is not currently its acceptance standard. This also excludes
mixed edits where Fortran answers are preserved and only JAX answers change.
Revisit this boundary when JAX is used in a host or a dedicated JAX BFB contract
is adopted; do not invent backend-specific marker variants in the meantime.

## Use the `-O0` preservation standard

CLUBB normally requires BFB at `-O0`; the current BFB regression tests use
that optimization level. A change that is BFB at `-O0` is considered BFB for
this policy even if it changes bits at `-O2` or another optimized level.
Compiler optimization can change floating-point results after small performance
refactors that preserve the mathematics and even the source order of operations.
Optimized-only differences do not by themselves justify the marker.

Note known optimized-only differences in the PR, with the compiler, optimization
flags, tested configuration and affected results, alongside the `-O0`
preservation evidence. Do not infer `-O0` BFB from mathematical equivalence or
unchanged source operation order; verify it with the relevant comparison.
Investigate unexpected differences and keep preservation claims scoped to the
actual evidence.

## Require a concrete numerical reason

Under the preservation standard above, use the marker only for an intentional
Fortran numerical answer change, with at least one of these justifications:

- **An affected test:** name the existing BFB test(s) that would fail without
  the declaration, the changed quantities and the relevant configuration.
  Distinguish an observed failure from a source-supported expectation.
  An explicit group such as "all Fortran BFB regression tests" is sufficient
  only when the effect actually spans that group. Check the configured cases,
  flags, compiler/backend and comparator; a job name alone is not evidence.
- **A precise condition outside current coverage:** identify the option values,
  code path and numerical calculation that changes answers. Explain which
  existing tests do not exercise it. For example, write
  `BIT_CHANGING only when OPTION=.true.` with the actual option name, and
  describe the affected host or standalone result. Include other necessary
  conditions such as PDF choice, input regime or backend. A green default run
  does not erase a known change in an untested branch.

If neither justification can be stated, investigate or report uncertainty;
do not add the marker merely because BFB has not been proved. Build/runtime
failures, infrastructure failures, console text and timers do not establish
a numerical answer change by themselves. Investigate unexpected failures
instead of declaring them expected.

Adding/removing a stat, reindexing a registry, or changing file metadata/schema
alone does not justify the marker. Existing common-variable BFB comparisons
can intentionally ignore stat membership changes. Check the actual comparator;
report these interface changes normally. Removing an option that no longer
has a supported run does not by itself show changed answers in a retained
configuration. A change to a retained stat's numerical values, or a diagnostic
calculation feeding back into the model, still needs the test-or-condition
decision above.

## Protect the remaining regression checks

Use the marker as sparingly as possible. It changes which answer differences
the BFB harness accepts and can conceal unrelated regressions in a large commit.
It is not a general waiver of builds, runtime checks or numerical validation;
some harnesses also flag declarations that produce no answer changes.

When a small answer-changing edit can be separated from a larger BFB-preserving
change, commit it separately. Keep the preserving edits unmarked and validate
them against their parent. Where practical, put preserving cleanup first and
the minimal numerical change afterward. Separate Fortran and JAX changes when
their preservation status differs and they can be reviewed independently.
In the PR, identify the marked commit and its affected tests/conditions.

Conditional wording documents scope for humans. The current
[BFB harness parser](https://github.com/larson-group/sys_admin/blob/master/python_nightly_test_suite/tests/bit_diff_test.py)
recognizes the marker in commit-message text without interpreting that condition;
it does not limit acceptance to a flag value or one named test. Use a space
after the marker for an explanation; the colon suffix is the harness's
retroactive commit-hash syntax, not a separator for prose.

## Check the final messages

Before publishing, review the actual commit messages and the planned merge or
squash message. Every declaration must have the test-or-condition justification
and the smallest practical commit scope. Quoting, backticks or negation do not
escape the parser: explanatory mentions of the literal marker in a preserving
commit message can still declare that commit. Keep such discussion in the
documentation or PR conversation instead.

Use [testing and reporting](../clubb-testing-and-reporting/SKILL.md) for evidence
and [PRs and commits](../clubb-documentation/change-messages/SKILL.md) for wording.
Keep verified BFB claims scoped to actual comparisons; a tolerance match or
numeric equality alone does not establish bitwise equality.
