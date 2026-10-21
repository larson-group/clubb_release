---
name: clubb-testing-and-reporting
description: Select relevant validation for a change and report Jenkins and local test evidence in PRs, commits and completion summaries.
---

# Choose validation and report what it establishes

Use when planning checks, finalizing a change, or writing its test results.
Test selection and reporting use the same question: what evidence exercises the
changed behavior and supports the claim? Keep technical acceptance requirements
in the relevant component/porting workflow; this skill does not replace them.
For the full sequence to get a branch ready to merge, use
[finalization](../clubb-finalization/SKILL.md). Standalone test selection or
reporting does not require that workflow.

For a reporting-only request, inspect recorded results rather than launching
new tests. For execution, stay within the task's authorized validation scope;
this skill does not require a full Jenkins suite for every change. Use
[Jenkins access](../clubb-jenkins/SKILL.md) for live jobs and
[pytest guidance](../clubb-pytests/SKILL.md) for focused suite ownership/admission.

## Choose checks that exercise the change

Read the diff, affected callers and existing test entry points. Identify the
routines, behavior, backends, cases and non-default options that could change.
Inspect pipeline commands and effective options, not just test names, to find
the relevant branch jobs. Prefer existing relevant Jenkins validation for an
integrated result; use focused local tests to investigate or cover gaps.

When numerical preservation is expected and a plausible concern, consider the
branch-versus-master comparison defined by
[run_bindiff_w_varying_flags_against_master/Jenkinsfile](../../../jenkins_tests/run_bindiff_w_varying_flags_against_master/Jenkinsfile).
Resolve its live branch job and inspect the compared revisions, flags, cases,
iteration limits and comparison criteria. Its pass covers that configuration,
not every compiler, backend or full-duration case. Follow
[numerical diagnosis](../clubb-numerical-diagnosis/SKILL.md) for attribution and
scientific comparison requirements. A tolerance pass does not establish
bit-for-bit (BFB) equality.
Do not require numerical comparison for changes that cannot affect results.

JAX correctness currently targets agreement with the corresponding Fortran
implementation, rather than preserving outputs from previous JAX versions.
Report JAX/Fortran comparison outcomes with their tested scope and limitations.
Use the canonical [bit-change policy](../clubb-bit-changing/SKILL.md) when
classifying answer changes and naming affected tests or untested conditions.
Explain numerical differences and unexpected failures regardless of whether a
declaration applies.

A focused test should exercise the changed routine or a meaningful edge/error
path. A large unrelated suite is not a substitute. For a regression fix, use
an active case that reaches the defect and check the baseline when practical;
do not claim the test catches the old bug without evidence.

## Decide when to wait

Running existing relevant Jenkins tests on the working branch can be part of
normal task validation without a separate user request for each run. Creating
or changing Jenkins tests requires a user request; follow
[Jenkins test authoring](../clubb-jenkins/test-authoring/SKILL.md) for that work.

During intermediate development, do not wait for every Jenkins test before
continuing independent work. Short checks can establish basic correctness while
larger comparisons take hours, especially when many jobs share resources.
Consider whether a result is needed for the next action and whether useful
independent work can proceed while it runs.

During intermediate development, start the relevant jobs, wait for fast checks
that meaningfully exercise the change, then continue work and revisit slower
results later. Retain queue/build links and tested revisions so pending checks
are not forgotten; later edits may require another run. Investigate failures
that undermine the next step's assumptions before relying on those assumptions.

For finalization, follow the
[finalization completion requirement](../clubb-finalization/SKILL.md#confirm-correctness-of-the-version-being-finalized):
wait for the applicable checks and review their results before declaring the
branch finalized. A fast pass does not replace required broader validation.
Any progress handoff or draft PR while jobs are still running must link them
as pending and state what remains unconfirmed; it is not completed finalization.
Never report launch or early stages as a completed pass.

## Verify the recorded evidence

For each result considered for the report, establish:

- The concrete job/build or local command, its result, and relevant stage(s).
- The source revision actually exercised; for branch/master comparisons, both
  revisions. A parameter naming a branch is not proof of its checked-out SHA.
- The relevant cases/options/backend and important exclusions, skips or limits.
- A permanent HTTPS URL to the specific Jenkins build, with a useful stage/log
  link when it helps explain the result. Do not use moving `lastBuild` links.

Confirm whether later changes affect the tested paths. A prior pass may remain
applicable after documentation-only edits; explain meaningful untested changes
rather than presenting an unrelated revision's green result as current evidence.
Distinguish success, failure, unstable, running, aborted, skipped and not run.
A queue entry or build aborted before the relevant stage is not a pass. An
overall failure can contain a successful relevant stage; state both accurately.
If results or links are unavailable, say what is known and what is missing;
never invent a build, outcome or URL.

## Prefer useful Jenkins evidence in a PR

Lead the Testing section with the most relevant Jenkins branch results that
were actually run, and link their specific builds with clear pass/fail status.
Highly relevant recorded passes must be referenced, especially a branch/master
comparison supporting an expected BFB claim or a job exercising changed code.
Leave out unrelated green jobs that do not add confidence in this change.

When Jenkins ran the same test, its result is the preferred PR evidence. A
local repetition usually need not be listed too. Include it if it tested the
change more rigorously or added coverage, such as a longer run, another backend
or flag, an exercised edge case or a stronger assertion; explain that addition.
This reporting preference does not prohibit useful local debugging or checks.
If no relevant Jenkins result exists, report the local evidence honestly and
state any important missing integrated validation.

Other local tests and pytests are useful when they add confidence beyond the
Jenkins results. Explain what was checked and why it matters. For example,
"The saturation-routine pytests passed, including its freezing-point boundary"
is more useful for that change than only "500 pytests passed." Counts can
supplement the explanation; do not use a large total in place of its purpose.
Only describe coverage confirmed by the actual tests/results.

When a Jenkins job or stage fails, check the same test on master with comparable
settings and record the tested revisions. Use an existing master result when
applicable; otherwise run the relevant comparison. A matching, understood master
failure need not be highlighted as a branch issue. Do not assume a failure is
pre-existing without evidence.

Report the relevant stages separately even when the overall job fails on master.
For a JAX change, for example, show the focused checks, standard comparisons and
other affected stages; a known failure in another stage does not replace their
results or make them irrelevant. Link the builds and distinguish stage results
from the overall status without repeatedly listing unrelated baseline failures.

Report every unexpected failure not established as a master baseline, including
its stage/command and effect on validation. Leave attribution open when evidence
is missing. A later successful rerun does not erase an unexplained failure.

## Adapt the report to where it is read

Use these evidence rules for PRs, commits and user-facing completion summaries.
Keep the most relevant result first, with a short explanation of what it
establishes and material gaps. PRs and summaries use clickable build-result
links. Commit messages remain concise but retain the key job/build, outcome
and result URL when Jenkins evidence supports the commit's claim. Do not carry
PR-wide validation onto an earlier commit that was not covered.

Local results should name the routine/feature and tested behavior, with the
command or test entry point when useful for reproduction. Identify the extra
rigor or missing Jenkins coverage that makes a supplemental check worth
mentioning. Preserve important failures and limitations while omitting redundant
success lists. A smoke pass, focused check and full comparison support different
claims; never shorten runs, weaken tolerances or hide failures to claim success.
