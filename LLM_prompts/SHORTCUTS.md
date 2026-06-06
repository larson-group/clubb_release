# LLM Prompt Shortcuts

Use this file as a repo-local index of reusable prompts.

When a user request appears to match one of these shortcuts:

1. Read the linked prompt.
2. Confirm that the prompt's goal and constraints match the user's current request.
3. If it matches, use it as task guidance.
4. If it only partially matches, use the relevant parts and state what does not apply.

## Shortcuts

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
