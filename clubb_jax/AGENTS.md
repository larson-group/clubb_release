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
[pytest workflow](../.agents/skills/clubb-pytests/SKILL.md). Focused checks live in
`pytests/`; new coverage goes in `pytests/auto_llm_generated_pytests/` pending
human review. JAX-specific real case runs/whole-driver workflows live in
this directory's `tests/`; shared comparisons remain in root `tests/`.
