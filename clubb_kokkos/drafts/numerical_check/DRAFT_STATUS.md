This is the complete supported source draft of
[`numerical_check.F90`](../../../src/CLUBB_core/numerical_check.F90).
It has not been compiled or run. Promotion and the maintained CPU Fortran versus
CPU Kokkos comparison belong to the integration owner.

All 12 Native procedures are present in Native order. The public
[`numerical_check.hpp`](numerical_check.hpp) preserves the eight Native entry
points. The existing `parameterization_check`, `check_negative`, `check_nan_1d`
and `check_nan_2d` definitions, plus the live public parameterization declaration,
are copied byte for byte from [`starting_live/`](starting_live/). The scalar
conservation residual is header-defined so device kernels can call it.

The full 1,889-line Native file and 466-line JAX file were read before drafting;
[`SOURCE_OUTLINE.json`](SOURCE_OUTLINE.json) records the contracts, branches,
comments and call order. Native remains the authority: this draft includes the
full PDF-parameter schema and settings checks that JAX reduced.

[`source_audit.json`](source_audit.json) records the separate source review:
94 newly ported calls and 371 arguments, 238 complete diagnostic payload values,
52 settings predicates, 11 parameter reductions and 17 device bodies. The audit
also checks enclosing IF branches for calls, writes, fatal assignments and
returns. These are source checks, not numerical evidence.
[`cleanup_audit.json`](cleanup_audit.json) records the later comment/format review.

See [`INTEGRATION.md`](INTEGRATION.md) for exact pending lower interfaces,
constants, compile gates and caller memory contracts, and
[`BFB_NOTES.md`](BFB_NOTES.md) for the numerical limits.
