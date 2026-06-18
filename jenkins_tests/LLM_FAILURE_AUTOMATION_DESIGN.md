# LLM Failure Analysis and Repair for Jenkins Tests

## What the current Jenkins setup implies

The Jenkins jobs under `jenkins_tests/` are mostly standalone declarative pipelines with:

- one or more `sh` steps that run compile/test commands
- an optional `cleanWs(...)` call in `post { always { ... } }`
- a `post { failure { ... } }` email hook using `emailext`

Two patterns matter for this design:

1. The jobs are not using a shared Jenkins library today, so any new behavior either needs to be added to each Jenkinsfile or routed through a common external script.
2. Some jobs clean the workspace in `always`, which means any failure-analysis hook that depends on workspace artifacts must run before cleanup or must copy artifacts to a persistent location first.

The BFB jobs are also special: they already delegate execution to an external wrapper:

- `/home/pub/jenkins_BFB_test_configs/sys_admin/python_nightly_test_suite/run_test.py`

That external wrapper is the best insertion point for richer per-step log capture for the BFB family.

## Recommended architecture

Use a two-layer design.

### Layer 1: deterministic artifact collection

Add a single post-test Python entry point, for example:

- `jenkins_tests/llm_failure_automation/llm_failure_report.py`

Its job is not to "reason" first. Its first job is to collect and normalize build evidence:

- Jenkins metadata: `JOB_NAME`, `BUILD_NUMBER`, `BUILD_URL`, `NODE_NAME`, `WORKSPACE`, `GIT_BRANCH`, `GIT_COMMIT`
- exit status / final Jenkins result
- stage logs captured during execution
- fallback console log if stage logs do not exist
- common output artifacts:
  - `output/**`
  - `build/**/CMakeFiles/CMakeError.log`
  - `build/**/CMakeFiles/CMakeOutput.log`
  - `*.log`, `*.out`, `*.err`
  - Python tracebacks and test summaries
  - compile output

This collector writes a small, stable bundle:

- `llm_artifacts/manifest.json`
- `llm_artifacts/evidence/*.log`
- `llm_artifacts/context.json`
- `llm_artifacts/report.md`
- `llm_artifacts/report.json`

### Layer 2: LLM diagnosis

After evidence collection, the script runs a diagnosis pipeline:

1. Local classifier:
   - identify failure type using regex and simple rules
   - extract the smallest useful snippets
   - rank candidate root-cause windows
2. Claude summarization:
   - summarize individual large logs if needed
3. Claude final diagnosis:
   - produce a report in strict JSON plus Markdown
   - cite evidence by file and line number
   - distinguish facts from hypotheses

This two-layer split is what makes the system autonomous in practice. The model should not receive raw unbounded logs as the first step.

## Why this is better than a single post-failure API call

A naive design would send the Jenkins console log directly to Claude after a failure. That will work sometimes, but it will fail operationally for large or noisy jobs:

- logs will exceed context limits
- unrelated warnings will dominate the prompt
- repeated failures will produce inconsistent diagnoses
- repair mode will not know which files or commands are actually relevant

The collector/classifier stage makes the output smaller, cheaper, and more reliable.

## Minimal Jenkins integration

### Short-term integration

For the first version, add a single post step to each Jenkinsfile and keep the existing email hook.

Use `unsuccessful` if you want to include unstable jobs, or `failure` if you want only hard failures.

Example:

```groovy
post {
  unsuccessful {
    script {
      sh '''
        set +e
        source /etc/profile.d/larson-group.sh || true
        python3 jenkins_tests/llm_failure_automation/llm_failure_report.py \
          --workspace "$WORKSPACE" \
          --job-name "$JOB_NAME" \
          --build-number "$BUILD_NUMBER" \
          --build-url "$BUILD_URL" \
          --node-name "$NODE_NAME" \
          --git-branch "$GIT_BRANCH" \
          --git-commit "$GIT_COMMIT" \
          --out-dir "$WORKSPACE/llm_artifacts"
      '''
      archiveArtifacts artifacts: 'llm_artifacts/**', allowEmptyArchive: true
    }
  }
  failure {
    script {
      emailext(
        to: 'messnermet@uwm.edu',
        subject: "${env.JOB_NAME} build ${env.BUILD_NUMBER} has Failed",
        attachLog: true,
        body: "${env.JOB_NAME} build ${env.BUILD_NUMBER} has failed. See the attached build log, the build results (${env.BUILD_URL}), and llm_artifacts/report.md."
      )
    }
  }
  cleanup {
    script {
      if ( "${env.JOB_NAME}" == "clubb_gfortran_build" ) {
        cleanWs(cleanWhenAborted: true, cleanWhenFailure: true, cleanWhenSuccess: true, cleanWhenUnstable: true)
      }
    }
  }
}
```

Two implementation notes:

- move `cleanWs(...)` to `cleanup` or another later hook so the report still has access to workspace artifacts
- the report step should be best-effort and should never hide the original build failure

### Better medium-term integration

Standardize shell execution through a wrapper so each stage produces a named log file.

For example:

- `jenkins_tests/llm_failure_automation/run_with_log_capture.sh`

Example usage:

```groovy
stage('Compile') {
  steps {
    sh '''
      jenkins_tests/llm_failure_automation/run_with_log_capture.sh \
        01_compile.log \
        "source /etc/profile.d/larson-group.sh && module load gcc netcdf-fortran && ./compile.py"
    '''
  }
}
```

The wrapper should:

- enable `set -o errexit -o nounset -o pipefail`
- create `llm_artifacts/evidence/`
- run the command with `2>&1 | tee ...`
- preserve the original exit code
- write a small metadata file containing stage name, command, exit code, start/end times

This is the most important change if you want good diagnoses.

## Suggested `llm_failure_report.py` behavior

The script should be a normal command-line program, not a Jenkins-specific Groovy script.

Suggested interface:

```text
python3 jenkins_tests/llm_failure_automation/llm_failure_report.py \
  --workspace "$WORKSPACE" \
  --job-name "$JOB_NAME" \
  --build-number "$BUILD_NUMBER" \
  --build-url "$BUILD_URL" \
  --node-name "$NODE_NAME" \
  --git-branch "$GIT_BRANCH" \
  --git-commit "$GIT_COMMIT" \
  --out-dir "$WORKSPACE/llm_artifacts" \
  --anthropic-api-key-env ANTHROPIC_API_KEY \
  --mode report
```

Suggested internal flow:

1. Read Jenkins metadata from args and environment.
2. Discover candidate evidence files.
3. If stage logs do not exist, optionally fetch the Jenkins console log or capture a local console surrogate.
4. Run a local failure classifier with patterns such as:
   - Fortran compiler error
   - linker error
   - Python traceback
   - unit-test assertion failure
   - SCM case diff failure
   - BFB mismatch
   - segmentation fault / abort / signal
   - missing module / environment failure
5. For each matched class, extract small windows around the highest-signal lines.
6. Build `context.json` with:
   - metadata
   - failure class
   - snippet list
   - likely failing command
   - candidate relevant files
7. Call Claude:
   - first pass only if logs are too large
   - final pass always uses the reduced context
8. Write:
   - `report.json`
   - `report.md`
9. Exit zero even if the API call fails, but write a local fallback report explaining why.

Suggested script shape:

```python
def main() -> int:
    args = parse_args()
    metadata = collect_metadata(args)
    evidence = collect_evidence(args.workspace, args.out_dir)
    classification = classify_failure(evidence)
    reduced_context = build_reduced_context(metadata, evidence, classification)

    if reduced_context.requires_chunking:
        chunk_summaries = summarize_large_artifacts(reduced_context)
    else:
        chunk_summaries = []

    try:
        diagnosis = run_claude_diagnosis(
            metadata=metadata,
            classification=classification,
            evidence=reduced_context,
            chunk_summaries=chunk_summaries,
        )
    except Exception as exc:
        diagnosis = build_fallback_report(metadata, classification, reduced_context, exc)

    write_manifest(args.out_dir, metadata, evidence, classification)
    write_json_report(args.out_dir, diagnosis)
    write_markdown_report(args.out_dir, diagnosis)
    return 0
```

Suggested module split:

- `llm_failure_report.py`: CLI entry point
- `collect.py`: artifact discovery and metadata capture
- `classify.py`: regex/rule-based failure classification
- `prompting.py`: Claude prompt construction and response validation
- `render.py`: Markdown and JSON report output
- `repair_handoff.py`: optional later bridge into autonomous repair mode

### Required output schema

The JSON report should be machine-readable. Example fields:

```json
{
  "job_name": "clubb_gfortran_build",
  "build_number": 12345,
  "result": "FAILURE",
  "failure_class": "fortran_compile_error",
  "confidence": 0.91,
  "summary": "Compilation failed in radiation_module.F90 after an argument list mismatch.",
  "root_cause": "A call site and callee signature diverged.",
  "evidence": [
    {
      "source": "llm_artifacts/evidence/01_compile.log",
      "line_start": 415,
      "line_end": 421,
      "reason": "Compiler diagnostic naming the file, routine, and mismatch"
    }
  ],
  "likely_files": [
    "src/Radiation/radiation_module.F90"
  ],
  "reproduction_steps": [
    "source /etc/profile.d/larson-group.sh",
    "module load gcc netcdf-fortran",
    "./compile.py"
  ],
  "suggested_next_actions": [
    "Compare the failing call site against the current subroutine signature",
    "Rerun the compile job after the signature fix"
  ],
  "notes": [
    "Evidence is from compiler output only; no runtime logs were produced"
  ]
}
```

This schema matters because it lets you do more than email a Markdown blob:

- send Slack or email summaries
- post a PR comment
- feed the result into a later repair agent
- track recurring failure classes over time

## Claude query pipeline

Use a strict pipeline rather than a single prompt.

### Step A: deterministic reduction

Done locally in Python:

- scan logs
- trim to high-signal windows
- remove repeated lines
- collapse long diff sections
- rank evidence by failure likelihood

### Step B: optional chunk summaries

Only if the reduced evidence is still large:

- summarize each large artifact independently
- preserve filename and line-number anchors
- save those summaries into `artifact_summaries.json`

### Step C: final diagnosis

Send Claude:

- build metadata
- failure class guess from the local classifier
- top evidence windows
- optional chunk summaries
- repository hints

Repository hints should be small and specific:

- common build commands
- known test families
- a short mapping from job family to likely code areas

Do not send the entire repository tree. Give only the parts relevant to the failing job.

### Prompt contract

Require the model to return:

- a factual summary
- a most-likely root cause
- alternative hypotheses if confidence is low
- evidence citations
- concrete reproduction commands
- concrete next debugging actions
- a statement of uncertainty when appropriate

Do not allow free-form speculation without citations.

## How to make the report path fully autonomous

"Autonomous" here should mean "no human intervention required to produce the report," not "no controls."

To make that operationally true:

1. Store the Anthropic API key in Jenkins credentials and inject it only into the report step.
2. Make the report step best-effort:
   - retries on network failure
   - bounded timeout
   - never changes the original build result
3. Archive the artifacts every time:
   - `report.md`
   - `report.json`
   - `context.json`
   - raw evidence snippets
4. Post the result automatically:
   - email attachment
   - Slack message
   - GitHub PR comment when the build belongs to a PR branch
5. Track recurring failures in a small datastore or even a flat file / database table:
   - failure class
   - affected job
   - commit SHA
   - whether the diagnosis was later confirmed

If you skip item 5, the system will still produce reports, but it will not improve over time.

## Recommended rollout

### Phase 1: report only

Add the report generator and artifact archiving. Do not change repair behavior.

This phase is highly feasible and low risk.

### Phase 2: report plus GitHub comment

If the Jenkins job is tied to a branch or PR, post `report.md` automatically as a comment.

This is still low risk.

### Phase 3: suggested patch only

Let Claude propose a patch, but do not apply it automatically. Save the patch as an artifact or PR comment.

This is a good intermediate checkpoint because it lets you measure usefulness before giving write access.

### Phase 4: autonomous repair branch and draft PR

Only after the report quality is stable should you allow automated code changes.

## Feasibility of autonomous fixing and PR creation

This is feasible, but only for a subset of failures.

### Good candidates

- broken Jenkinsfile logic
- missing or wrong script paths
- simple Python harness regressions
- deterministic compiler errors
- straightforward CMake/configuration issues
- small API drift between call sites and definitions

### Medium candidates

- deterministic runtime crashes with clear stack traces
- missing environment/module setup on the Jenkins node
- output-format mismatches in scripts

### Poor candidates

- scientific correctness regressions
- BFB/numerical mismatch failures
- tolerance-related output drift
- non-deterministic parallel failures
- failures caused by external machine state or flaky infrastructure

For CLUBB specifically, the "poor candidate" list matters. Many failures in numerical model code are not safe to repair from logs alone.

## How repair mode should work

Do not implement repair mode as "one Claude prompt that edits the repo."

Use an agent loop on a dedicated repair worker:

1. Detect a qualifying failed build.
2. Clone the failing commit into a fresh worktree.
3. Reproduce the failure with the exact Jenkins command.
4. Run the report pipeline first.
5. Start a repair agent with:
   - repository access
   - the structured report
   - permission to edit files
   - permission to run a limited command allowlist
6. Require the agent to:
   - reproduce the failure before patching
   - make the smallest plausible fix
   - rerun the failing test
   - run one or two related smoke tests
7. If validation passes:
   - create a branch
   - commit with a bot identity
   - open a draft PR
8. Attach:
   - the original failure report
   - the validation commands
   - the post-fix result

### Required safety controls for repair mode

- use a dedicated bot GitHub account
- only create branches, never push to protected branches
- open draft PRs only
- enforce per-run time and token budgets
- restrict shell commands to an allowlist
- require reproducible failure before patching
- require successful rerun before PR creation
- store every prompt, patch, and validation log as artifacts

Without these controls, the system will create noisy or unsafe PRs.

## Best insertion points by job family

### Regular Jenkins jobs in this repo

Add:

- the post-failure report step in each Jenkinsfile
- the log-capture wrapper around stage commands over time

### BFB jobs

Add:

- top-level post-failure report call in the Jenkinsfile
- richer command-level logging inside the external `run_test.py`

This keeps the Jenkinsfile simple while giving the BFB suite better evidence.

## Recommended first implementation

If the goal is to get value quickly, implement exactly this first:

1. Create `llm_failure_report.py`.
2. Add `post { unsuccessful { ... } }` report generation to 2-3 high-value jobs first:
   - `clubb_gfortran_build`
   - `clubb_python_test`
   - `clubb_BFB_e3sm_flags_gfortran_test`
3. Move any `cleanWs(...)` that would remove evidence before the report runs.
4. Archive `llm_artifacts/**`.
5. Keep the existing email notification.

Only after that should you standardize the stage log wrapper and consider repair mode.

## Bottom line

Report generation with Claude is highly feasible now and can be made fully autonomous with modest Jenkins changes.

Autonomous repair with branch creation and draft PRs is also feasible, but it should be treated as a second project with stricter controls. It will work well for infrastructure, script, and deterministic compile failures. It will be much less reliable for numerical-model regressions and BFB-style failures unless you require strong validation before a PR is opened.
