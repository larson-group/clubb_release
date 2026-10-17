---
name: clubb-jenkins
description: Run and monitor Jenkins tests on a branch, inspect failures or hangs, and access live jobs through the HTTPS API.
---

# Run and inspect Jenkins tests

Use when asked to run tests, inspect builds or failures, or work with live jobs.
For creating or changing test definitions, use
[Jenkins test authoring](test-authoring/SKILL.md); its requirements also
cover coordinated live-job changes. The
[Jenkins README](../../../jenkins_tests/README.md) explains setup and reproduction.
Read the current `jenkins_tests/` definitions and inspect the live job inventory,
configuration and recent build logs through the Jenkins HTTPS API. The login
script sets up each user's credentials; this skill guides agents through job
operations without browser clicking or a server SSH login.

## Login and API access

Use `https://carson.math.uwm.edu/jenkins/` directly over HTTPS. Set up the intended
Jenkins account's username and API token with the standard-library-only
[login script](../../../jenkins_tests/jenkins_api_login.py). No SSH access is required.

Before job operations, look for `~/.config/jenkins/api-credentials` on the machine
running the API requests. It contains one `USERNAME:API_TOKEN` line. Read it
inside Python without printing it, require owner-only file permissions on POSIX,
and verify the identity with `GET whoAmI/api/json`. Require an authenticated
identity whose name is not `anonymous`; an HTTP 200 alone is insufficient.
Existing credentials in this format work without running setup again.

If the file is missing, malformed, or authentication is rejected, tell the user
to run this from the repo root in their own interactive terminal:

```bash
python3 jenkins_tests/jenkins_api_login.py
```

They generate a named API token in the intended account's Security/Configure
page, or use a token assigned to them for that account, then enter its username
and token at the prompts. Rerun setup to switch accounts or tokens.
The token prompt is hidden. The script verifies the login before atomically
saving or replacing the plaintext credential file with mode `600` on POSIX.
Never ask the user to paste a token into chat or put it in a command argument.
Resume after setup and recheck identity. For broad file permissions, ask the
user to run
`chmod 600 ~/.config/jenkins/api-credentials` before using the file.

For API calls, use Python's `urllib.request` with an HTTP Basic Authorization
header constructed from the saved line. The login module's `JENKINS_URL`,
`CREDENTIALS_FILE`, and `verify_credentials()` can be imported from
`jenkins_tests.jenkins_api_login` when running from the repo root. Keep tokens and
authentication headers out of command text, logs and artifacts. Use HTTPS and
reject redirects for authenticated requests. Jenkins may return HTTP job/queue
links: rebuild those links under the configured HTTPS base after checking the
host and `/jenkins/` path before following them. Do not forward credentials to
another origin. API-token requests do not need browser cookies or CSRF crumbs.

Read jobs/configuration with `api/json` and `job/NAME/config.xml`; parameterized
builds use `POST job/NAME/buildWithParameters`. Read the live job's parameter
definitions before supplying its branch parameter (normally `BRANCH`). Encode
job names and parameter values, including nested `job/FOLDER/job/NAME` paths.
Follow the returned queue item until a build is assigned, then monitor its
`api/json` and `consoleText`. Inspect the checked-out revision and report the
HTTPS build link. A request to test unpushed local changes needs a reachable
commit before Jenkins can check it out; do not claim local edits were tested.

If login succeeds but an operation returns 403, check that user's Jenkins
permissions rather than asking for SSH access or another person's credential.
Report DNS, TLS and connection errors as connectivity problems; do not replace
valid credentials or bypass HTTPS verification to hide them. API operations are
limited by the same account permissions as the browser. Use the relevant API
endpoints for authorized job/view changes.

Jenkins references: [authentication](https://www.jenkins.io/doc/book/system-administration/authenticating-scripted-clients/)
and [remote API](https://www.jenkins.io/doc/book/using/remote-access-api/).

## Account permissions

A token has its owning account's permissions; this workflow can use restricted
or administrator accounts. A successful identity check does not establish what
that account can do. Relevant permissions are:

- `Overall/Read` and `Job/Read` for ordinary job/result access; `Job/Build` to
  launch builds and `Job/Cancel` to stop them.
- `Job/ExtendedRead` for read-only job configuration access. It may need enabling
  before an administrator can grant it; it does not grant configuration edits.
- `Job/Create` to create a job with its complete initial configuration through
  the API; `Job/Configure` to change an existing job's saved configuration.
  Creation alone does not allow later corrections or edits.
- `Job/Delete` for job deletion; administration and other operations have
  separate permissions. See [Jenkins permissions](https://www.jenkins.io/doc/book/security/access-control/permissions/).

Before promising live job changes, establish the needed access for the affected
job/folder using available role information or permission-gated reads. If access
is unverified, qualify the plan; do not change live objects merely to probe it.
If blocked, identify the denied operation and needed permission, continue useful
repository work, and report the site action as pending. Readable configuration
or permission to build does not prove permission to edit. Account permissions
and the user's authorization for a particular action are separate requirements.

## Scope and inventory

Distinguish branch, production, release, BFB and host-model jobs. The user's
request determines which can change or run. A branch validation request does
not authorize production configuration changes. Existing relevant branch jobs
can run as normal task validation without a separate request for each build;
creating or modifying Jenkins tests requires the user's request under
[Jenkins test authoring](test-authoring/SKILL.md). If asked to run a Jenkins job,
trigger that actual job and retain its queue/build URL and configured revision;
a manual shell approximation does not establish that Jenkins works.

Cancel only builds or queue items launched for the current task unless the user
has authorized stopping others. Track their queue/build links; a shared account's
username does not establish which person or task started a build.

Include disabled jobs in the inventory; disabled does not mean obsolete.
Check the configured source revision and pipeline path rather than inferring
coverage from the job name. For naming, coverage changes and coordinating site
updates, use [Jenkins test authoring](test-authoring/SKILL.md).

## Verification

Use [testing and reporting](../clubb-testing-and-reporting/SKILL.md) when choosing
jobs and reporting their results. Retain permanent links to relevant builds,
verify their tested revisions/coverage, and report unexpected failures. Follow
its [guidance on when to wait](../clubb-testing-and-reporting/SKILL.md#decide-when-to-wait):
use fast relevant results to continue intermediate work, revisit slower jobs
later, and consider waiting longer when final validation needs their results.

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
