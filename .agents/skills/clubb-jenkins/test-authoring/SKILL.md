---
name: clubb-jenkins-tests
description: Create, edit, consolidate, rename or remove CLUBB Jenkins tests and pipelines, with reproducible commands and coordinated live-job changes.
---

# Create and maintain Jenkins tests

Use for changes to test definitions and coverage in
[jenkins_tests/](../../../../jenkins_tests/README.md). Read the relevant existing
Jenkinsfiles and script entry points before choosing the change. For live job
inventory, API access or validation builds, also use
[Jenkins access](../SKILL.md). A request only to run or diagnose a
job does not need this authoring workflow.

## Change tests only at the user's request

Create or modify Jenkins test definitions or their configured coverage only
when the user requests that work. A feature, option or bug fix alone does not
authorize adding, removing, renaming or changing Jenkins tests. Suggest a new or
updated test when it would usefully cover the change, especially for a new
feature or option, and wait for the user's go-ahead before authoring it. Authorization already given
in the conversation is sufficient; do not ask again.

Running existing relevant Jenkins tests on the working branch is appropriate
as part of normal validation and does not need a separate request for each run.
Use [testing and reporting](../../clubb-testing-and-reporting/SKILL.md) to choose
checks and decide when their results are needed. Live job configuration changes
still require coordination as described below.

## Keep steps easy to reproduce

Prefer a short sequence that a user can copy from the Jenkinsfile into a
terminal at the repository root: select the environment, call a script with
explicit options, or run a simple Linux command. Keep non-trivial calculations,
case selection, parsing and test logic in their existing script owners rather
than adding Groovy or embedded Python/shell programs to the pipeline. Reuse an
existing owner before introducing another wrapper.

Prefer literal options and case names over environment variables or Groovy
variables that hide the command. Avoid new environment variables when a script
option or direct path suffices. Jenkins-provided metadata, credentials,
resource controls and unavoidable worker setup may still require variables.
Explain what these exceptions supply and why they are needed; do not remove
working configuration solely to avoid a variable.

Follow the [shared CLI conventions](../../../../run_scripts/README.md#shared-command-line-conventions)
and [script workflow](../../clubb-script-cli/SKILL.md): single-dash descriptive
options, explicit cases on run lines and the existing dependency bootstrap in
script entry points. Avoid redundant defaults. Let Pyplotgen select cases from
its input folders unless the test specifically needs a filter. Indent command
continuations and use concise comments for unusual commands. Prefer inline
`-override` JSON when it makes case settings visible.

## Explain the test and unusual steps

Start each new or changed Jenkinsfile with a `//` comment block explaining:

- What behavior is tested, which variants/cases matter, and what the comparison
  or pass/fail result establishes.
- The compiler, environment, inputs and working directory needed to reproduce it.
- A simplified set of real terminal commands for reproducing the relevant check.

Keep reproduction commands aligned with the scripts and options actually used.
For plotting, external BFB suites or other workflows that cannot be reduced to
a useful short example, explain the limitation and point to the actual entry
point and required inputs rather than inventing a simplified equivalent. Use
[file descriptions](../../clubb-documentation/file-descriptions/SKILL.md) for the explanation's
structure and [code comments](../../clubb-documentation/code-comments/SKILL.md) for local rationale;
Jenkinsfiles use Groovy `//` comments.

Comment unusual environment variables, GPU locks, worker restrictions,
machine-specific deployment and other unavoidable orchestration where they
occur. Explain the reason, relevant resource/host and how it affects reproduction.
The reader should understand the workflow from the comments and commands.

## Parallel work and failure behavior

Run independent tests or backend/compiler variants in parallel where practical.
Give each parallel branch its own checkout, build and output files; do not race
on a shared build directory, generated input or comparison output. Respect
actual dependencies and serialize access to a shared GPU or deployment target.
Explain a required lock or host label rather than forcing parallel execution.

Keep GPU case workers at one; use the agreed resource limit for long CPU tests
and retain Pyplotgen's own automatic process policy. Reuse the configured
SCM wipe/reclone lifecycle rather than adding `cleanWs` to pipelines.

Preserve independent failure reporting when combining checks, and keep manual
aborts effective. Do not catch an intentional abort as an ordinary failure and
continue expensive work. Do not add retries or weaken cases/tolerances to mask
an underlying failure. Retain appropriate results/artifacts and notifications.

## Coverage, names and Jenkins site changes

Add backend coverage to the job that owns the native/Fortran equivalent; driver
lifecycle coverage belongs in `clubb_driver`, for example. Distinguish production,
branch, release, BFB and host-model jobs. Preserve intended case/compiler/flag
coverage when combining or removing tests, including disabled jobs. Standard
full-run coverage does not imply coverage of unmaintained or multi-year cases.
Branch names use `clubb_branch_*`; do not append `_branch` or broadly rename
unrelated job families. Release jobs follow their configured release revision.

Whenever adding, removing or renaming a test, remind the user that the
corresponding jobs/configuration on the Jenkins site must also change. A repo
commit does not create, rename or remove live jobs. List the required site
changes, including affected production/branch jobs, SCM Jenkinsfile paths,
parameters, schedules and views where relevant. Before promising to apply them,
check [account permissions](../SKILL.md#account-permissions); being able to author
a Jenkinsfile does not establish access to create or edit its live job. Keep
blocked site actions explicit even when repository work can proceed.
For removals, map the old job to its replacement and preserve history/backups
as appropriate.

Coordinate these changes with the user before applying them. Editing a
Jenkinsfile does not by itself authorize live job/view changes. Running existing
relevant branch jobs follows the normal validation scope above.
Once the user has authorized specific site actions, use the HTTPS API through
[Jenkins access](../SKILL.md) and report what was actually changed;
do not ask again for actions already covered by that authorization.

## Check and report the change

Use [testing and reporting](../../clubb-testing-and-reporting/SKILL.md) to select
validation for the affected behavior and report the relevant build links,
supplemental local coverage and unexpected failures.

Check the real script options, copied reproduction commands, pipeline syntax
and affected coverage. Use the live Declarative validator when access permits;
it validates syntax, not the test's scientific coverage. Run only the builds
within the user's requested scope. Record the exact job/build and source
revision, stages reached, and whether validation was a smoke check or a full
pass. Report pending site changes separately from completed repository work.
