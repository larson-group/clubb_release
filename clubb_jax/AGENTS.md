# JAX port work

For Fortran-derived source edits or audits in this directory, read
[`LLM_prompts/jax_porting_workflow.md`](../LLM_prompts/jax_porting_workflow.md)
and its faithful-port standards. The structural completion audit applies even
when numerical tests pass: preserve routine/file order, definitions and calls,
comments/dividers, readable spacing, diagnostics and non-default branches.

Use `README.md` for the current runtime/support boundary. The supported
standalone is JAX-owned; do not introduce Fortran/F2PY fallbacks. Scope an
explicitly approved numerical-core rewrite separately from its public interface.
