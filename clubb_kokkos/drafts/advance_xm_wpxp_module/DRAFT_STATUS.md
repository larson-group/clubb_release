# Mean/flux advance draft status

The full supported `advance_xm_wpxp_module.F90` is translated in the native-named
C++ source/header. The public 107-argument parent and all 17 private procedures
are present in native order. The source body performs its own assembly,
coefficient calculations, orchestration, updates, diagnostics and stats; only
public lower dependencies remain explicit f2k calls.

The supported scope follows the actual JAX exclusion of sponge damping. The
three native sponge boundaries are marked in the source and require upstream
rejection of sponge-enabled cases before integration. Full native diagnostics
are retained, including the 82-argument fatal field dump; JAX's reduced dumps
were not copied.

This is an isolated source draft, not an integrated or compiled port. No model,
unit, oracle or timing test was run. No BFB or performance result is established.
No production source, CMake, generator, test or root worklog was edited.

Completed separate source audit and cleanup:

- Read the complete 6472-line native module, actual JAX scope, existing grid,
  matrix, stats and array owners, and public lower routine contracts.
- Compared all 18 signatures/order and both genuine optionals. Checked all 151
  retained source-order subroutine calls and 1744 non-solver arguments, including
  sections, indices, field-specific tolerances and named optional mappings.
- Compared 328 explicit assignment operation trees independently, preserving
  source arithmetic grouping. Re-read the 12 MIN/MAX assignments, whole-array
  sections, solver permutation and documented mask-definedness choices.
- Checked 192 kernel bodies for owning-UDT capture, host reads/I/O, nested
  launches and bridge calls. Device kernels use component Views/scalars only.
- Reviewed loop independence, per-column clipping scan order, all stats/error
  branches, all 161 non-solver diagnostic write/label sequences and ownership.
  The two solver-error writes are reviewed at their explicit adapter. Ran the formatter;
  delimiter, whitespace and 100-column checks pass.

`source_audit.json` and `formula_audit.json` record the source-only evidence.
Their scripts never execute physics. `lower_contracts.json` records actual
native public dependencies; the adjacent declaration-only f2k header does not
implement those bridges. `prepare.py` recreates an initial translation and must
not overwrite any numerical fixes made after handoff.

Root owns compilation, incremental live hookup and the maintained model runner.
The regular numerical gate is CPU Fortran versus CPU Kokkos with zero tolerance;
GPU comparisons are occasional milestones under the latest user policy. All
options, scalar/column configurations, error paths, native-definedness edges and
GPU behavior remain runtime untested for this draft.
