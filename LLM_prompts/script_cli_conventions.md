# CLUBB script and CLI conventions

Use when adding/changing CLUBB script arguments or their callers, including
Jenkins commands. This records the script-argument cleanup conventions; it does
not require renaming unrelated APIs, third-party flags or structured MCP fields.
The option contract is owned by [run_scripts/README.md](../run_scripts/README.md#shared-command-line-conventions);
check that section and the relevant script help instead of maintaining a second
option inventory here. [utilities/README.md](../utilities/README.md) documents
shared namelist behavior; [setup_python_venv.py](../utilities/setup_python_venv.py)
owns the bootstrap. [tests/README.md](../tests/README.md) documents comparison
harness controls.

## Public option names and meanings

| Concern | Convention |
| --- | --- |
| CLUBB options | One dash and descriptive snake_case, e.g. `-output_dir`. |
| Independent process/case concurrency | `-workers`; default half the available logical CPUs, at least one. Reuse `tuner.system_defaults.default_max_workers`. |
| Destination paths | `-output_dir`, `-output_file`, `-output_root` as appropriate, rather than `-out_dir` or vague `-out`. |
| Parameter inputs | `-params_file`, `-param_ranges`, or the existing specific name (`-silhs_params_file`, `-flag_config_file`). |
| Multi-column generation | `-multicol NUM_OR_SPEC` through the canonical namelist utilities. |
| Model/tuner batch width | `-batch_size`; this is model columns/candidates, not worker-process concurrency. |
| Other parallelism | `-threads` is OpenMP threads; `-process_counts` is a benchmark sweep. Neither is an alias for workers. |

Tuner public destinations use `-output_job_dir` / `-output_run_dir`; internal
`-job_dir` selects an existing job. Preserve these documented meanings rather
than applying a blind path-name replacement.

Use the existing parser and `-help` to verify the actual accepted surface before
changing examples/callers. Preserve intentional positional arguments and
third-party command options such as `nvidia-smi --format`. Do not mechanically
replace every `--` in a shell command. New compatibility aliases need a concrete
consumer/migration reason; do not reintroduce obsolete spellings as a precaution.

CPU worker defaults apply to comparable worker controls, not every integer
called jobs/batch size. GPU run stages use one case/process worker to avoid
device contention; this does not prohibit vectorized device batches. Pyplotgen
keeps its own automatic max-process policy; do not add process-control options
or impose the CPU half-count convention on it. Explicit user/site limits win.
Avoid redundant defaults in Jenkins commands, but retain explicit serialization
when it explains a GPU/resource constraint.

## One owner for shared behavior

`utilities/create_case_namelist.py` owns namelist options and interpretation.
Reuse `add_namelist_arguments()` and `parse_forwarded_args()`; wrapper scripts own
only their own orchestration options and forward namelist/run options to the
existing owner. The same applies to JAX selection/runtime setup: SCM wrappers
forward `-jax[=VALUE]` to the JAX launcher instead of duplicating backend/venv
logic. Executable scripts that need third-party Python dependencies should call
the existing `ensure_python_venv` bootstrap before those imports. Preserve
arguments/environment across relaunch. Standard-library-only paths need no
forced bootstrap, and importing a library module must not restart its caller.
JAX environments remain owned by the JAX launcher.

`-override` accepts a simple all-case override, inline JSON, or a JSON file.
JSON supports `all` defaults followed by case-specific values; the case wins
on a duplicate key. Parse this once in the namelist utility, not in each runner.
Keep execution inputs, quoting and forwardable argument values intact. Reject
invalid/unowned arguments clearly rather than silently dropping them.

Prefer the existing owner module for small helpers. A new helper module is
justified by distinct shared behavior, not a few thin forwarding wrappers.
Do not grow defensive retry/compatibility layers to hide a broken owner contract.
A small lexical adapter can be necessary to keep a bare `-jax` from consuming
a positional case; preserve that explained grammar without duplicating runtime
policy. Model bounds/flag changes belong in the canonical Fortran validator and
its existing Python mirror, documented in `utilities/README.md`, not a second
parameter policy in a runner or physics kernel. Local numerical guards may
still be needed.

## Changing and checking a CLI

Search repo callers (Jenkinsfiles, wrappers, tests, docs) and known external
consumers. Update the full authorized chain. Report required host-model/BFB
migration separately; do not alter live configs or host repositories beyond the
user's requested scope. Check forwarding, inline/file overrides, malformed
inputs, and affected default behavior with focused parser/runtime checks.
Compare any failing validation with the baseline before calling it a regression.
Describe options in plain language: what they control, accepted values, default,
and any meaningful restriction. Indent shell line continuations relative to the
command, and add concise one-line comments for unusual commands or constraints.
For longer scripts, show meaningful stages and a final case/settings/result
summary; keep detailed diagnostics in logs and provide a usable failure follow-up
command when available. Direct comparison callers should consume the existing
compact structured status, not parse human log wording. Higher-level diagnostic
tests may inspect logs for the specific evidence they exercise.
