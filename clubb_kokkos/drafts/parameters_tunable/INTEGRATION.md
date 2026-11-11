# Tunable parameter behavior

This draft implements all eight agreed supported procedures in Native order:
`set_default_parameters`, `check_parameters_api`,
`get_parameter_hard_bounds_api`, private `check_parameter_hard_bounds`,
`calc_derived_params`, `init_clubb_params_api`, private `pack_parameters` and
private `unpack_parameters`. The two Native tuner readers are outside the actual
JAX/standalone scope and remain explicitly omitted. The complete Native file,
actual JAX mirror and driver contracts were read before the source outline.

The existing `nu_vertical_res_dep` mirror remains the only schema owner.
Parameter matrices, derived diffusivities, `lmin` and `mixt_frac_max_mag` remain
owned by `ClubbDriver` or its caller. Native has no mutable module instance here.
The padded 105-name catalog and Fortran-flat 210-value bounds are immutable
`constexpr std::array` objects; there is no global Kokkos allocation or cache.
The existing `parameter_indices.hpp` already contains zero-based array indices.
Physical grid metadata remains Native one-based and is translated only when
accessing the height Views.

Merge into the existing generated header with the established preserved-region
mechanism. Keep the authoritative type region intact. The procedure order is:

1. `procedure_includes.hpp.inc` in the header include section.
2. Existing generated `nu_vertical_res_dep` type.
3. `constants.hpp.inc`, `private_declarations.hpp.inc`, then
   `procedure_definitions.hpp.inc`, inside `namespace clubb_kokkos`.
4. Add permanent `parameters_tunable.cpp` to the core library. Its top-level
   includes are in `cpp_includes.inc`; put `init_clubb_params_api.cpp.inc` inside
   `namespace clubb_kokkos`.
5. Merge `clubb_api_module.hpp.inc` at the Native single/multi-column derived
   wrapper positions, with `api_includes.hpp.inc`.

The initializer's same-named HostSpace rank-2 overload is out of line because
`driver_io.hpp` includes generated UDT conversion headers, which in turn include
`parameters_tunable.hpp`. Keeping the transient reader in the `.cpp` avoids this
include cycle. A thin same-named generic output overload preserves a host,
strided or device caller's allocation using the shared array owner. The actual
contiguous HostSpace driver uses the fixed overload without a device round trip.
This is a language/I/O boundary, not a second initialization implementation.
The other templates and three source-private helpers retain their Native body
order in the header; private forward declarations satisfy their earlier callers.

Default values are written by a Kokkos column kernel using the output Views'
execution space. Initialization owns one call-local HostSpace matrix and 105
named contiguous column subviews; packing uses an independent call-local matrix.
The private pack kernel selects the output matrix's execution space. These
column arguments must have homogeneous compatible View types/accessibility;
actual initialization satisfies that contract. The private unpack routine stages
one parameter vector once, then writes 105 scalar references. No global copy
slots, per-scalar Fortran calls or new permanent data type are needed.

Extend only the existing `generate_driver_io.py` owner with the staged
`generate_driver_io.py.inc` insertion, before its individual-control section.
It derives all 105 member declarations and their exact order from Native
`init_clubb_params_api`'s `/clubb_params_nl/`. `driver_io.hpp.inc` and
`driver_io.F90.inc` are generated source previews for review, not new generators
or separate compilation units. The required reader contract is
`f2k::read_clubb_params_nl(runfile, ngrdcol, iunit, <105 initialized arrays>, read_status)`.
It performs only file open/namelist read/close with call-local values and the
existing wire descriptor owner. No Native default, pack, unpack or whole-init
routine is called. Required `iunit` remains present; the existing I/O owner's
`NEWUNIT` convention replaces Native `OPEN(unit=iunit)` at this exact boundary.
An empty or all-blank filename skips reading, matching Native padded character
comparison. A nonzero read status raises a C++ fatal exception after the service returns;
compiler-specific Native I/O failure text and list-directed spacing are unverified.

Replace reached whole-procedure bridge calls together:

- Driver initialization currently calls `f2k::init_clubb_params_api` once for the
  full parameter-set matrix; replace it with the permanent API.
- Both driver derived-value calls use `f2k::calc_derived_params_api`; replace them
  with the new permanent overload.
- `numerical_check.cpp::check_clubb_settings_api` currently calls
  `f2k::check_parameters_api`; use the permanent check and its existing header.
- Retire those three obsolete bridge entries when their final callers move.
  Private pack/unpack/check helpers are never exported as Native bridge APIs.
- The driver still obtains `params_list` through its existing environment
  snapshot. Its printing/stats metadata can use the immutable C++ catalog when
  that snapshot member is retired; preserve all 32 padding bytes where the
  existing descriptor contract expects them. This is metadata, not mutable state.

Validation preserves all-column traversal and per-column fatal updates. Hard
bounds run before the four C6 equalities; existing error codes are not reset and
no fatal return suppresses later diagnostics. NaNs retain Native comparison
semantics. Derived values allocate all seven nu arrays only when `nu1` is absent,
write `lmin`, and stop on any invalid minimum before mixture or diffusivity math.
Prescribed/type-3/type-1/type-2 spacing branches and all seven nu mappings remain
in source order. Invalid `grid_type` without prescribed spacing leaves Native
`avg_deltaz` undefined; the existing setup-grid gate excludes that direct-API
boundary. No JAX fallback or invented initialization is added.

The single-column wrapper represents Native sequence association with an owned
one-column parameter matrix. Native creates such a local matrix but passes its
original vector to the explicit rank-2 dummy; values and ordering are identical.
Integer powers 2 and 3 are expanded mechanically. No FMA, reassociation barrier,
new clipping or speculative BFB arithmetic workaround is introduced.

`source_audit.py` records exact signatures, 315 default/pack/unpack assignments,
24 derived assignment expression trees, four equality predicates, all catalog
entries, both wrapper contracts and I/O membership/call order. `CLEANUP_AUDIT.json`
is a separate final grouping/comment/whitespace pass with exact pre/post-format
executable-token comparisons. These are owner source checks. No independent
review, build, model execution, CPU/GPU comparison or BFB claim is made; all
runtime branches await compatible-group integration.
