# Jenkins tests

This directory contains the pipelines Jenkins uses to build CLUBB, run tests
and publish results. Each test subdirectory holds a `Jenkinsfile`. This README
explains how to set up API access, find a test's commands and coordinate changes
to its live Jenkins job.

## Find and run a test

Open [CLUBB Jenkins](https://carson.math.uwm.edu/jenkins/) to select a job and
inspect its stages, console output and results. Production jobs generally use
the name of their pipeline directory. Corresponding branch jobs use the
`clubb_branch_*` prefix and normally accept a `BRANCH` parameter. Check the live
job's configuration and parameters before running it; release and external
host-model jobs have their own revision/coverage rules.

For example, [clubb_python/Jenkinsfile](clubb_python/Jenkinsfile) runs the
Python API, shared unit checks and Python/Fortran comparisons.
[clubb_basic_full_run/Jenkinsfile](clubb_basic_full_run/Jenkinsfile) runs the
standard cases with two compilers in separate parallel branches. See the
[test harness README](../human_reviewed_tests/README.md) for the underlying test commands and
[run-script README](../run_scripts/README.md) for ordinary case execution.

## Set up an agent's HTTPS access

[jenkins_api_login.py](jenkins_api_login.py) saves the credentials an AI coding
agent needs to run and inspect jobs using your Jenkins account. You do not need
a server SSH login or an open browser session for subsequent API requests.
Run setup in your own interactive terminal, from the repository root:

```bash
python3 jenkins_tests/jenkins_api_login.py
```

Generate a named API token in the intended account's Security/Configure page,
or use a token assigned to you for that account. Enter that account's username
and token at the script's prompts; token input is hidden. The script verifies
your identity before saving one `USERNAME:API_TOKEN` line in
`~/.config/jenkins/api-credentials`, outside the repository. This is plaintext storage with owner-only mode `600` on POSIX.
Rerun setup to switch accounts or replace a token; a failed login preserves
the previous file. Existing credentials in this format continue to work.

For setup help without prompting:

```bash
python3 jenkins_tests/jenkins_api_login.py -help
```

The script uses only Python's standard library. It sets up access; it does not
submit tests. The [Jenkins access skill](../.agents/skills/clubb-jenkins/SKILL.md)
guides an agent through selecting jobs, supplying the branch, following the
queue/build and reading results. API requests use the saved token with your
account's Jenkins permissions; signing in does not grant additional permissions.

Reading and building normally use `Overall/Read`, `Job/Read` and `Job/Build`;
stopping builds requires `Job/Cancel` on the relevant jobs. `Job/ExtendedRead` adds
read-only configuration access. Creating jobs and editing existing jobs require
separate permissions, so an agent should check access before promising site
changes. See the skill's [account-permission guidance](../.agents/skills/clubb-jenkins/SKILL.md#account-permissions).

## Reproduce a check in a terminal

Read the pipeline's opening description and `sh` steps for its environment,
compiler, case inputs and actual commands. Opening descriptions on new or
changed tests should include copyable reproduction commands. Older pipelines
may still require reading their steps. Run commands from the repository root
unless the pipeline specifies another directory.

For example, the commands in
[clubb_budget_balance/Jenkinsfile](clubb_budget_balance/Jenkinsfile) are:

```bash
source /etc/profile.d/larson-group.sh
module load gcc netcdf-fortran
./compile.py
python3 human_reviewed_tests/check_budget_balance.py .
```

This example uses the UWM workers' environment setup. Elsewhere, prepare the
corresponding compiler and dependencies using the
[build instructions](../BUILD.md).

A plotting comparison or BFB (bit-for-bit) job can require additional reference
files, hosts or setup. Follow its actual pipeline and entry point rather than
assuming a short command list recreates it. For example, the
[BFB pipeline](clubb_BFB_gfortran_test/Jenkinsfile) uses the shared external
`/home/pub/jenkins_BFB_test_configs/sys_admin/python_nightly_test_suite/run_test.py`
harness and its named configuration. Reading logs through HTTPS does not make
unpublished files on a worker available as build artifacts.

## Create or change a Jenkins test

Use the [Jenkins test-authoring skill](../.agents/skills/clubb-jenkins/test-authoring/SKILL.md)
for the full requirements. Keep pipeline steps simple enough to copy into a
terminal: environment commands, script calls with explicit options and simple
Linux commands. Put substantial test logic in its script owner. Prefer direct
arguments to environment variables, run independent work in parallel with
separate files/workspaces, and explain necessary variables, locks and
machine-specific behavior. Describe what is tested and how to reproduce it at
the top of the Jenkinsfile.

Repository definitions and live Jenkins jobs are separate. Adding, removing or
renaming a test also requires coordinated changes on the Jenkins site: its
jobs, SCM pipeline paths, branch parameters and views as relevant. Include that
reminder and the exact required actions when proposing a change. An agent can
apply authorized site changes through HTTPS, but must coordinate with the user
before doing so; editing the repo is not permission to change live jobs.

## Choose checks and report results

The [testing-and-reporting skill](../.agents/skills/clubb-testing-and-reporting/SKILL.md)
explains how to select checks that exercise a change and report their evidence.
In a PR, prefer links to relevant Jenkins branch builds and their actual
outcomes. Mention local repetitions when they add rigor or coverage, and explain
what extra focused tests checked. Report unexpected failures and important
untested behavior rather than presenting only a large passing-test count.
