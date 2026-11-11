# Remapping draft

[remapping_module.cpp](remapping_module.cpp) ports all 18 procedures from
[the native owner](../../../src/CLUBB_core/remapping_module.F90), in native order.
The public entry points are `remap_iv_from_var_name` and the 13-argument
`remap_vals_to_target`. The latter returns an owned `(column,target cell)`
`LayoutLeft` View. Its generic overload retains the source field's memory space
and accepts const, strided inputs through the existing array owner.

The draft has owner source/cleanup review plus Serial and CUDA Kokkos compilation. It has
no whole-model numerical evidence: all reached remapping callers are currently
inside retained Fortran owners. Root explicitly accepted that boundary rather
than introducing a reverse callback or porting higher state solely to test it.

## Promotion

1. Copy the two native-named `.cpp`/`.hpp` files into `src/CLUBB_core/`.
2. Apply only [integration.patch](integration.patch): one CMake source entry and
   two selected native constants in the existing type generator. Regenerate the
   selected model flags with that owner; preserve all other current selections,
   source-preserved type regions and lifecycle changes.
3. Remove draft-status wording from the header when promoted. Keep the exact
   adaptations/limitations below. No procedure bridge removal or physics bridge
   addition is needed: this draft has no current C++ call site to replace.
4. Label the module source-audited/compiled until a real caller reaches it. Use
   the maintained whole-model CPU comparison when that outer call is ported.

`cons_ullrich_remap=1` and `ppm_remap=2` come from
[model_flags.F90](../../../src/CLUBB_core/model_flags.F90), lines 94–97. Existing
physical/rational constants, the grid schema, scalar interpolation API and
host debug query supply everything else. No private native symbol is declared
or called in science code.

## Callers retained in Fortran

- `atex_tndcy` remaps total-water and temperature forcing when
  `l_add_dycore_grid` is true; its enclosing forcing call remains native.
- `apply_time_dependent_forcings_from_dycore` uses native-held dimension/forcing
  objects and remains inside the native time-dependent-input owner.
- The stats output writer remaps ZT/ZM fields for `l_different_output_grid` with
  NetCDF output enabled. Stats remain the approved native backend.
- Adaptive-grid remapping remains outside the supported driver scope.

The first three are fixed-grid uses too; remapping is not classified as dead
merely because adaptive grids are excluded. The exact inventory is recorded in
[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json).

## Language and execution adaptations

- Native grids omit `INTENT`; the complete public body only reads them. C++ uses
  `const grid&`. The private optional operator is a pointer defaulting to null;
  every public input stays required. Supplying it still computes masses and
  ordered pressure recurrences before the operator selection/debug checks.
- Fields are column first. Private PPM coefficients adapt native
  `(coefficient,column,level)` to `(column,level,coefficient)`; the four entries
  remain AA, AL, AR, A6. Profile subviews map native `i1:i2` onto local columns.
  Latitude storage maps the native lower bound `jfirst-ng_s` to zero. `k0` keeps
  its native level-zero slot explicitly.
- A scalar-cell `kmppm` overload is the device/coefficient-shape adaptation;
  the original four-argument array form calls it. Independent interior limiter
  calls are fused into one cell kernel. It never reads unwritten A6 before the
  native `kord=4` Lin branch initializes it.
- Singleton edge/field dimensions replace native scalar sequence association in
  the PPM wrapper and consistency diagnostic. All allocations remain call-owned.
  The public PPM wrapper still calls one profile at a time and reuses its target
  scratch between columns. Matrix-vector and PPM full-cell sums remain ascending
  folds. Native `SUM` in the integral diagnostic is represented by an ascending
  fold; exact compiler association has not been compared.
- `map1_ppm` searches independent columns together, preserving each column's
  ascending moving search. Host code selects the first failing native column
  and reports it before any output interpolation. No host bridge or diagnostic
  executes inside a device kernel.
- Five global `zero`/`one_half` operands passed to device `Kokkos::min/max`
  are copied into same-kind `real_t` prvalues because that API binds scalar
  references. Values, operators and association are unchanged; this is device
  portability, not a BFB arithmetic workaround.
- Native error stops become C++ exceptions with the same messages; writes use
  host stderr with all source data. Native list-directed spacing is not claimed.
  Debug2 still checks the same data in the same order. No numerical threshold,
  safe divisor, search clamp, FMA or source-definedness repair was introduced.

## Source scope and definedness

Both Ullrich and the public PPM `kord=4` path are implemented. `steepz` is
unconditionally called before the Huynh gate and is part of the active path.
The draft additionally retains the small private Huynh `kord>=7` branch, without
adding a public order setting or claiming numerical coverage. The native
`km>32` wind sponge is retained even though JAX omits it. Native sequential mass
integration is used rather than JAX's cumulative-sum formulation.

The original input assumptions remain: increasing, covering density splines,
compatible grids, at least four source PPM cells and bracketed monotone pressure
edges. With only three cells, the native provisional AL loop is empty while
the top stencil reads AL(3); no new input guard or output repair was added. A missing PPM bracket prints its original payload and returns before
writing `q2`. The diagnostic's implied `kk=1..km+1` loop can read beyond `pe2`
when `kn<km`; this source-definedness issue is retained and documented beside
that loop. The public reused target scratch can consequently be unspecified on
its first failed profile or retain the preceding output on a later failure.
No replacement output is invented.

Copied source/target count labels and density/mass units were corrected only in
comments where the implementation establishes their meanings. All 18 procedure
names, 118 native arguments, 44 selector names and 14 explicit source call sites
are retained. See [source_audit.json](source_audit.json),
[formula_audit.json](formula_audit.json) and [cleanup_audit.json](cleanup_audit.json).
These artifacts record owner reviews. A distinct agent completed a full-body
review in `/tmp/clubb_remapping_independent_source_audit.json`; its three
documentation/format/device-boundary findings are corrected, and the final-hash
recheck is recorded there. No independent runtime claim is made.
