---
name: clubb-finalization
description: Prepare a CLUBB branch for merge by auditing its changes, verifying current test evidence, resolving readiness decisions with the user, and committing, pushing and creating or updating a draft PR.
---

# Finalize a branch for merge

Use when the user asks to finalize a branch or get it fully ready to merge.
A request only to review code, run a check or draft a message does not invoke
this whole workflow. Finalization includes the commit, push and draft-PR handoff
below, subject to any narrower user instructions; it does not include merging.

"Finalize right away", "go ahead and finalize" and similar requests mean to
proceed through the whole workflow without routine review pauses. They do not
waive auditing or testing. Substantive, unaccepted issues still trigger the
stopping rule below; urgency does not authorize silently fixing or deferring them.

## Recover intent when useful

Optionally begin with available session history and the local worklog. Compare
user requests and accepted plans with the actual branch to catch accidentally
unfinished work. Distinguish those obligations from abandoned approaches,
superseded requests and suggestions the user never approved. Include meaningful
omissions in the audit; do not silently expand the task to revive old plans.

## Audit the whole branch before runtime validation

Establish the intended PR base and current branch state, including staged and
unstaged changes and intended new files. Review the complete branch diff against
its base, not only the last commit. Preserve unrelated user work. Check the
current target branch for conflicts or changes that affect the audit; resolve
routine conflicts within scope and revisit affected checks after integration.

Read the changed implementation and the helpers, callers and interfaces needed
to understand it. Select applicable standards through [AGENTS.md](../../../AGENTS.md),
including faithful porting, coding/formatting, CLI and host compatibility where
relevant. Audit logic, boundary/error behavior, implementation choices, coverage
and integration. Look for dead code, unused imports, obsolete paths and temporary
workarounds introduced or made unnecessary by this branch. Passing tests do not
replace a porting or implementation audit.

Use [documentation guidance](../clubb-documentation/SKILL.md) for comments,
headers, READMEs and change messages. Check factual claims against implementation
or recorded evidence, including examples, defaults, options, numerical/BFB and
performance claims. Check the existing PR description and relevant agent-authored
comments for stale claims; keep the current description accurate and clarify
material corrections when posting is authorized.

Fix minor wording, formatting and unambiguous cleanup within scope, then review
the resulting diff. Apply the stopping rule below as soon as a substantive issue
is found; do not launch validation already known to be doomed.

## Confirm correctness of the version being finalized

Testing must finish before finalization is complete. Use
[testing and reporting](../clubb-testing-and-reporting/SKILL.md) to select and run
applicable checks, including existing relevant Jenkins branch jobs, then wait for
and review their results. Local smoke tests do not replace those jobs. Follow
component acceptance requirements; a full unrelated suite is not mandatory. Use
[Jenkins access](../clubb-jenkins/SKILL.md) for live jobs and permission checks.
Finalization does not authorize creating or changing Jenkins tests or live jobs.

Only genuinely trivial prose changes, such as README wording or comment typos
that cannot change behavior, may omit runtime/Jenkins tests; review the diff and
state why the exception applies. Documentation changes to commands, defaults or
behavioral contracts do not qualify merely because they are in Markdown or
comments. If no existing Jenkins job exercises a change, run the applicable
owner checks and explain that coverage limit rather than running unrelated jobs.

Confirm that Jenkins results relied on exercised the latest branch version
being finalized, with the intended comparison baseline, cases and options.
Verify actual checked-out revisions rather than trusting a branch parameter or
an earlier green build. Run missing checks or rerun stale ones. If subsequent
edits affect the tested behavior, repeat the affected audit and validation.

Jenkins needs a reachable commit to test new changes. If necessary, finish the
audit and available local checks, then commit/push a validation checkpoint before
launching Jenkins; that checkpoint does not mean finalization is complete.
Continue independent work while slow checks run, then return to review the
results before completing finalization. A checkpoint or draft PR with tests
still running is progress, not completed finalization. Required results that
are pending, unavailable or inconclusive prevent claiming readiness; report
the gap rather than treating launch or a smoke pass as a full pass.

## Stop for substantive, unaccepted problems

If the audit or validation reveals a clear unintended issue that the user has
not already discussed and accepted, pause finalization before further substantive
implementation changes or the final PR handoff. Examples include incorrect logic,
a meaningful porting/compatibility defect, a missing requested behavior, a
materially false correctness claim, or a failure that would break a relevant
Jenkins test. Gather enough evidence to explain the issue; do not silently fix
it, weaken checks or redefine success to finish.

Tell the user what was found, the supporting evidence, its effect on merge
readiness and the proposed remedy. Before substantive edits, verify that the
user accepted that remedy or explicitly authorized working through issues;
permission to investigate alone does not satisfy this checkpoint. Diagnose and
explain the issue first, then wait for acceptance of a fix, deferral or limitation.
Repeat affected checks after the accepted change.

This checkpoint applies to unexpected, non-trivial failures during late testing
as well as finalization. Earlier implementation may allow routine problem-solving
within its authorized scope. A previously accepted limitation need not be
approved again, but keep it visible. Trivial corrections such as a missed rename,
a broken README extension, a comment typo or small formatting cleanup may proceed;
mention useful cleanup without turning every nit into a stop.

## Commit, push and hand off a draft PR

After the audit and required validation have completed and been reviewed, and
substantive failures have been resolved or explicitly accepted under the stopping
rule, review the final diff and working tree. Apply the
[bit-change policy](../clubb-bit-changing/SKILL.md)
to commit separation and audit the actual commit/merge messages before publishing.
Commit intended work and push the branch. If no changes remain
to commit, retain the existing commits; do not manufacture an empty commit.
Verify the remote head and that the test evidence still applies to that head.

Find the existing PR for this branch before creating another. Create new agent
PRs as drafts and keep agent PRs in draft state; do not mark ready or merge as
part of finalization. Use [PRs and commits](../clubb-documentation/change-messages/SKILL.md)
for wording, evidence and comment attribution.

For an existing PR, update its title/description when the final scope, behavior
or evidence has changed. A concise progress comment can help explain significant
additions since review when posting that update is authorized; do not rely on
comments alone when the description has become misleading. Do not post a comment
merely to announce that finalization ran. Comments posted through a user's
account end with the agent signature specified in the writing guidance.

Verify the PR's base, head, draft state, final diff and description after the
update. Report its link, relevant test outcomes, accepted limitations, minor
cleanup worth mentioning and any remaining human/site actions. Distinguish
ready for human merge review from pending validation or unresolved work.
