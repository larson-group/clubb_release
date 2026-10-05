# Agent Instructions

Before substantial work, check [LLM_prompts/SHORTCUTS.md](LLM_prompts/SHORTCUTS.md)
and read only the workflow that fits. Determine scope from the current request
and conversation; ask when an unresolved choice changes the outcome. Shortcut
selection does not itself authorize source edits, external writes, goals, or
broader validation. Explicit user instructions take priority.

## Repository conventions

- Fortran is the authority for faithful ports. Passing tests does not replace
  the routine/branch/signature/comment/formatting audit in the porting workflow.
- For script/CLI work, use
  [script_cli_conventions.md](LLM_prompts/script_cli_conventions.md): single-dash
  descriptive options, `output_*`, and `-workers` for process concurrency with
  a half-logical-CPU default. Keep worker counts distinct from batch width and
  preserve external-tool flags and documented GPU/plotting exceptions.
- Follow [pytest_workflow.md](LLM_prompts/pytest_workflow.md): fast specific
  checks belong in `pytests/`; new agent coverage goes in
  `auto_llm_generated_pytests/` pending human review. Real cases/full workflows
  belong in their component's `tests/`; reserve root `tests/` for shared
  workflows and keep Jenkins additions within the requested scope.
- Reuse the existing utility/parser/state owner. Keep UI/runner adapters thin;
  justify a new abstraction by actual shared behavior. Remove obsolete paths
  in the authorized change instead of adding silent fallback/retry machinery.
- Write user docs around purpose, a usable quick start, then needed detail.
  Update relevant sections of the existing document; keep implementation detail
  in its owner and define scientific symbols/units when they matter.
- Compare failures with the starting revision before attributing them to a
  change. Report supported, tested and untested behavior accurately.
- Preserve unrelated work and machine-local artifacts. Inspect the diff before
  staging; do not include ignored worklogs, outputs or local symlinks in commits.

## Session Worklog

Use repository-root `worklog.md` as the durable record for substantial work.
Read recent entries before continuing an investigation. Append a dated milestone
with the objective, decisions/assumptions, files, verification/results and remaining work;
never rewrite earlier history. Keep it concise and exclude secrets and large
raw output. This file is ignored: updating it does not mean committing it.

For an approved durable convention discovered during work, update its narrow
canonical guidance as part of the authorized change, or record an explicit
follow-up if documentation changes are outside scope. Keep task-specific choices
in the worklog; do not promote a question or one-off exception to a repo rule.
