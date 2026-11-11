This draft implements all three Native hydrometeor/precipitation lifecycle
routines, preserving their source inventory and order. The full 191-line Native
and 136-line actual JAX files were read before the outline. The existing generated
scalar and collection types remain the authoritative schema.

Merge `procedure_includes.hpp.inc` at top level and `procedure_declarations.hpp.inc`
inside `clubb_kokkos` after the generated types in the existing same-named header.
Compile this folder's `.cpp` as the procedure owner. Root must add preserved
procedure regions to the existing generator without duplicating the mirror types.

`init_hydromet_pdf_params` resets all 18 fields in Native order. The 14 fixed-size
Native array components use the existing C++ View fields: allocate missing storage
once, preserve matching allocations, and reject malformed already allocated
shapes. The existing collection owner supplies the canonical maximum of eight
hydrometeors. A thin same-named collection overload calls that owner's `zero()`;
it represents Native's repeated per-cell initialization and introduces no second
field inventory. The current C++ driver already owns and allocates this collection.

`init_precip_fracs_api` implements IntentOut by releasing the caller's previous
fields, allocating three `(column,thermodynamic-level)` Views in source order,
setting metadata and calling the reset. `zero_precip_fracs_api` preserves those
allocations and writes all three fields in Native per-cell order with column-first
LayoutLeft access. The copied Native argument comment incorrectly described
precipitation fractions as hydrometeor parameters; target comments use their
actual dimensionless fraction meaning.

Replace matching driver `f2k::init_precip_fracs_api` and
`f2k::zero_precip_fracs_api` calls with these routines. The existing driver
collection `.zero()` may use the same-named initializer overload to mirror the
Native call; its allocation owner and layout stay unchanged. Retire the two
corresponding bridge manifest entries after their final callers are replaced.
There is no new lower scientific bridge, constant, type or global state.

`SOURCE_AUDIT.json` covers the 18 reset fields, 14 fixed storage contracts,
three precipitation allocations/metadata and the reset kernel. Separate
`CLEANUP_AUDIT.json` records Native argument groups, complete descriptions,
formatting and exact executable-token preservation, with unused includes removed
separately. No FMA, association or other deliberate BFB adaptation was added.
Independent review, compilation and all initialization/reset branches remain
pending compatible-group CPU zero-tolerance validation.
