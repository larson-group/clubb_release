---
name: clubb-porting-jax
description: Use JAX arrays, functional state, compiled control flow, batching and autodiff without carrying inefficient source-language loops or state into a port.
---

# JAX representation and transformations

Use when JAX is source or target with the
[porting strategy](../../SKILL.md). These rules describe reusable language/runtime
considerations. Current model support, launch commands and project-generated
documentation belong to the implementation's maintenance owner.

## Decide state, types and shapes early

Use explicit array dtypes and precision configuration. Choose immutable model
objects/pytrees with dynamic array leaves and genuinely static metadata. Return
updated state rather than mutating hidden objects. Source-like public arguments
can coexist with an explicit context/state adaptation; preserve their mapping.

Python state read during tracing can be frozen into the compiled function.
Do not use mutable globals for dynamic model fields or parameters. Shape and
static-argument changes may specialize; runtime-value changes intended to reuse
compilation should actually do so. Measure compilation reuse rather than infer
it from a decorator. Avoid variable-length Python containers in compiled work.

Keep computational math in JAX. NumPy conversion, file I/O and opaque backend
calls belong at explicit host boundaries. Repeated `np.asarray`, `.item()`,
`float`, `int` or `bool` extraction of traced values can break transformations.
Use the established shared boundary owner rather than local converter helpers.
Read [pytrees](https://docs.jax.dev/en/latest/pytrees.html) and
[JIT](https://docs.jax.dev/en/latest/jit-compilation.html) for the installed API.

## Translate work, not the source loop syntax

- Use array expressions or `vmap` over independent columns/samples. Choose
  logical axes deliberately; source storage order is not automatically the
  appropriate batch/scan axis. Avoid Python loops performing one `.at` update
  per level, which can produce huge unrolled compiled programs.
- Preserve true vertical/time/solver recurrences using `lax.scan` or appropriate
  structured control flow. Do not apply `vmap` to dependent iterations.
- Keep static configuration branches in Python where appropriate. For traced
  conditions use `lax.cond`, masks or other structured primitives; do not branch
  through `bool(jnp.any(...))` or escape data with Python early returns.
- Whole-slice `.at` updates for bands, boundaries or fixed assembly are useful;
  avoid a blanket prohibition on updates. Keep source temporaries and operation
  order where they affect diagnostics or the numerical contract.
- Make invalid branch inputs safe before arithmetic. `jnp.where` does not make
  an invalid expression in its other operand harmless; NaNs/infinities can also
  contaminate gradients. For example, substitute a safe positive operand before
  evaluating a logarithm and then apply the intended output mask.

See [automatic vectorization](https://docs.jax.dev/en/latest/automatic-vectorization.html),
[scan](https://docs.jax.dev/en/latest/_autosummary/jax.lax.scan.html) and
[control flow](https://docs.jax.dev/en/latest/control-flow.html).

## If gradients are part of the goal

Check meaningful derivatives incrementally, alongside primal model comparisons.
A finite return or a constant zero gradient can hide a disconnected path.
Use controlled inputs and appropriate directional/finite-difference checks when
needed. A small probe does not establish all-physics or long-run differentiability.

Dynamic early termination can block reverse-mode differentiation. A static-bound
`scan`/`fori_loop` with a carried active/converged mask may be appropriate:

```python
def step(carry, _):
    state, active = carry
    candidate, converged = safe_iteration(state)
    state = select_active_leaves(active, candidate, state)
    return (state, active & ~converged), None
```

This is pseudocode for an iteration that tests convergence after an update;
match the source's exact test/update order. Freeze converged elements rather
than continuing their numerical updates. Preserve final errors/status, counters,
statistics and random draws; safe operands are needed even on masked work.
Only apply this design when the iteration bound and derivative meaning are
justified. It can cost more execution but enable the requested feature.
For dynamic bounds, `fori_loop` may lower to `while_loop`, whose ordinary reverse
mode is unsupported; see [fori_loop](https://docs.jax.dev/en/latest/_autosummary/jax.lax.fori_loop.html).

Use the normal runtime path for feature probes. A legacy bridge may be useful
while migrating, but cannot establish differentiation through the delegated
operation. Keep host effects outside pure compiled science or represent them
through explicit diagnostics/status; inspect callback transformation limits.

## Numerical and performance checks

Test representative compiled paths early even if throughput is not the goal:
avoiding eager dispatch and massive Python unrolling can make the port's own
test cycle much faster. Do not mechanically decorate every helper; structured
primitives already compile their bodies. Check startup/compile and steady-state
cost separately, synchronizing device work at measurement boundaries.

XLA may fuse/reassociate/simplify operations, changing rounding even within JAX.
Do not assume `jit` preserves eager bits or that every language pair can/cannot
be BFB. Follow the agreed [numerical contract](../../validation/SKILL.md) and
investigate the first difference before tolerating accumulated drift. Preserve
algorithm/permutation parity separately from random-stream identity.
See [numerical FAQ](https://docs.jax.dev/en/latest/faq.html#jit-changes-the-exact-numerics-of-outputs)
and [benchmarking](https://docs.jax.dev/en/latest/benchmarking.html).
