# JAX port work

For Fortran-derived source edits or audits in this directory, read
[`.agents/skills/clubb-jax-port/SKILL.md`](../.agents/skills/clubb-jax-port/SKILL.md)
and its faithful-port standards. The structural completion audit applies even
when numerical tests pass: preserve routine/file order, definitions and calls,
comments/dividers, readable spacing, diagnostics and non-default branches.

Use `README.md` for the current runtime/support boundary. The supported
standalone is JAX-owned; do not introduce Fortran/F2PY fallbacks. Scope an
explicitly approved numerical-core rewrite separately from its public interface.

For pytest additions, maintenance and CI coverage, read the
[pytest workflow](../.agents/skills/clubb-pytests/SKILL.md). New agent tests go in
`auto_llm_generated_tests/`: focused checks in `pytests/unit_tests/`, checks
against Fortran through the Python API in `pytests/f2py_unit_tests/`, and
JAX-specific real case runs/whole-driver workflows in `integration_tests/`.
`human_reviewed_tests/` holds tests a human has reviewed; shared comparisons
remain in the root `human_reviewed_tests/`.
