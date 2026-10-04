# Change or validate Jenkins tests

Use for Jenkins job/pipeline maintenance or a request to run named Jenkins tests.
Read the current `jenkins_tests/` definitions and inspect the live job inventory,
configuration and recent build logs through the available authenticated tools.
Use current connection discovery; do not hard-code a historical SSH tunnel,
port, token or account into this workflow. Keep credentials out of artifacts.

## Scope and inventory

Distinguish branch, production, release, BFB and host-model jobs. The user's
request determines which can change or run. A branch validation request does
not authorize production configuration changes. If asked to run a Jenkins job,
trigger that actual job and retain its queue/build URL and configured revision;
a manual shell approximation does not establish that Jenkins works.

Include disabled jobs in the inventory. Disabled does not mean obsolete.
For removal or consolidation, map repo definitions to live jobs and replacements;
remove superseded jobs only when deletion is in scope. Preserve the intended
case/compiler/flag coverage. Branch jobs use the `clubb_branch_*` prefix,
with production names mapped to their repo definition; do not append `_branch`
or broadly rename other families during an unrelated repair. Release jobs follow
the release revision, so a master merge alone does not migrate them.
Standard full-run coverage and an expanded list of
unmaintained or multi-year cases are different validation scopes.

## Pipeline changes

- Reuse the configured wipe/reclone lifecycle; do not add `cleanWs` to pipelines.
- Independent parallel stages may acquire their own node/executor, checkout,
  build and output. Confirm isolation rather than sharing mutable workspaces.
- Use current CLI conventions and the shared dependency bootstrap in the script
  entry point. Avoid copying venv setup into each caller to repair one script.
- Keep GPU run stages at one worker. Use the approved job/resource limit for
  long CPU tests; Pyplotgen retains its own automatic process policy.
- Indent continuations and explain unusual commands with concise one-line
  comments. Prefer inline `-override` JSON when it makes case settings visible.
- Keep intentional aborts effective; do not catch them as ordinary failures and
  silently continue expensive work. When combining tests, independent failures
  should remain visible without erasing which stage failed.

## Verification

Compare the changed behavior with the starting revision, including already
failing jobs. Record the source revision, exact job/build and stages reached.
An unchanged name or a running queue item is not proof of model coverage. For
an agreed smoke audit, confirm startup and relevant stages; distinguish that
from a completed full-duration pass. Abort/rerun only within the user's scope.

Before declaring a hang, inspect fresh console output, queue/executor state,
process progress and stage timing. Correct the actual code/pipeline cause when
it is in scope. Do not introduce retries or weaken cases/tolerances to hide it.
Report required live configuration or external host/BFB migrations separately
from repository changes, and apply them only when authorized.
