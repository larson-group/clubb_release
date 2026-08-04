# LLM Prompt Shortcuts

Use this file as a repo-local index of reusable prompts.

When a user request appears to match one of these shortcuts:

1. Read the linked prompt.
2. Confirm that the prompt's goal and constraints match the user's current request.
3. If it matches, use it as task guidance.
4. If it only partially matches, use the relevant parts and state what does not apply.

## Shortcuts

### Update Host Models After CLUBB Changes

Use when a CLUBB change affects a host-consumed Fortran/C API or interface: an exported `_api` routine's signature or semantics, public types/constants, generated wrappers, or host-facing flags/configuration. This workflow is for host-owned compatibility changes only. Never modify, copy, synchronize, reformat, or otherwise include vendored CLUBB/SILHS source in a host-model PR; source synchronization is a separate process. Do not use this shortcut for Dash/MCP changes or Python-only wrapper changes unless the host model directly consumes the changed interface. If host repositories are not provided, audit the branch and report the exact host-owned call-site change needed; do not clone or edit external repositories implicitly.

Prompt:

- `LLM_prompts/update_host_models_after_clubb_changes.md`

### Fix Python API

Likely use when the user asks for things like:

- fix the Python API
- update the Python driver after Fortran changes
- make Python/JAX drivers match Fortran
- repair f2py wrappers after refactors
- get `run_python_vs_fortran_cases.py` passing
- get `run_jax_vs_fortran_cases.py` passing

Prompt:

- `LLM_prompts/update_python_api_and_drivers.md`

Before using, confirm whether the user wants:

- Python API only
- Python driver only
- JAX driver too
- full comparison suite passing
- just compile/tests/smoke tests

### Port Underlying Fortran

Likely use when the user asks for things like:

- port underlying Fortran to another language
- re-port a stale Python, JAX, or other file from Fortran
- make a ported file match the underlying Fortran source file
- similarize a target-language port against the Fortran source
- remove target-only helpers, aliases, optionals, or reordered logic
- make routine calls, comments, or argument lists match the Fortran file

Prompt:

- `LLM_prompts/port_underlying_fortran_to_other_languages.md`

Before using, confirm whether the user wants:

- investigation only, or source edits now
- one target file only, or all language mirrors
- exact structural matching, or an idiomatic target-language rewrite
- validation only, focused tests, or full comparison suites

### Format Fortran Routines

Likely use when the user asks for things like:

- format a new Fortran subroutine
- clean up a routine interface
- extract a helper routine in CLUBB core
- make call-site formatting match a routine definition
- apply CLUBB routine formatting rules
- document routine arguments and local variables

Prompt:

- `LLM_prompts/fortran_routine_formatting.md`

Before using, confirm whether the user wants:

- formatting-only changes
- behavior-preserving refactor plus formatting
- full call-site and wrapper updates
- source code edits now, or only a review/checklist

### Work Through CLUBB Dash

Likely use when the user asks to:

- connect to, attach to, use, open, or control the Dash dashboard
- compile, run a case, inspect output, make or save plots, navigate a tab, or
  configure/launch/inspect a tuning job
- investigate a result, especially when the dashboard can show the relevant
  profiles, contours, run console, Tune settings, or activity log
- create, update, or view a static investigation report

Prompt:

- `LLM_prompts/dash_app_workflow.md`

Before using, determine whether there is one unambiguous open dashboard with a
live browser view.  Do not silently start another dashboard when none is open,
and do not guess between multiple instances.  Read the prompt before taking
dashboard actions; it defines the required connection, persistence, report,
and fallback behavior.
