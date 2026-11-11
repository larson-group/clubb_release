# Positive-definite flux-adjustment draft

This draft implements the complete Native `pos_definite_adj` routine, with its
10 required arguments and four caller-owned changed arrays. JAX supports this
same algorithm. `SOURCE_OUTLINE.json` predates implementation; source-expression
and separate cleanup records accompany the draft. The files are outside the
build and have no compilation or numerical evidence.

Promote `pos_definite_module.cpp` and `.hpp` to the core directory and add the
source to the core target before replacing the matching call in the reviewed
mean/flux draft. No new public lower contract or selected constant is needed:

- Existing C++ `ddzm(nzm,nzt,ngrdcol,gr,Real2D)` returns the derivative on
  thermodynamic levels. The Native `flux_lim-flux_np1` array expression is
  materialized once at that existing API boundary.
- Existing host-only `f2k::clubb_at_least_debug_level_api(int)` queries Native
  debug state, followed by the original `Correcting flux` stdout message.
  Its call-site TODO can disappear when the error-code owner is ported. There
  is no Native call inside a device kernel.
- Existing `grid`, `array_types`, `eps`, `zero_threshold` and `zero` remain
  the owners of layout, staging, metadata and values.

All physical and scratch fields use `(column,level)` LayoutLeft. Concrete
execution-space overloads consume the core grid; generic overloads use the
existing shared array owner for host or strided callers and restore all four
outputs to their original allocations. Native scratch is per call. Explicit
value copies preserve the old field and final flux; a View alias would change
those semantics. No global storage or new state accessor is introduced.

The limiter uses the exact Native physical loop bounds and signed grid step,
flattened by ordinal. It does not substitute a hard-coded ascending interior.
Local Native indices remain one-based until their View accesses. Its lower
limiter term deliberately uses `dz_over_dt(k_zt)`, even when the field and
outward flux use `kp1`. Boundary copies use storage endpoints `1,nzm`.
The two column ANY checks are computed separately, before and after the field
update, and gate entire-column flux/field budgets respectively. Both false
branches write zero. No floating reduction or vertical recurrence is reordered.

Full Native grid cartoons, scientific citation, equation sequence and the
untested-negative-at-nz FIXME remain. The second copied budget comment now names
`field_pd`, and the scratch description calls `flux_lim` the limited flux,
matching the actual assignment. Same-kind constant values are materialized for
Kokkos min/max device references; their values and argument order are unchanged.
No speculative FMA, reassociation or signed-zero/BFB workaround is present.

Compilation, CPU/GPU execution, forward/reversed grids, clipping/non-clipping
columns, debug output and strided callers are all untested. Parent integration
owns the maintained CPU Fortran-vs-Kokkos zero-tolerance model gate and occasional
GPU comparisons; a source review does not establish those results.
