# Hydrometeor metadata and correlation setup

This source-only draft implements all seven Native procedures in definition
order: two private mapping/default helpers, name lookup, initialization,
correlation setup, symmetry assertion and matrix display. JAX supports the
same algorithmic scope; the draft retains Native error objects and full
diagnostics. Compilation and numerical validation are pending.

The existing generated `corr_varnce_module.hpp::hm_metadata_type` remains the
only metadata type. Its 26 fields and current converters are untouched.
The current generator replaces that whole header, so promotion must establish
a preserved generated-type region before adding permanent declarations, as
already done for `grid_class`. Merge the provided fragments in this order:

1. Required includes: `array_types.hpp`, `err_info_type_module.hpp`, `<array>`
   and `<optional>`, in addition to the existing generated includes.
2. `parameter_types.hpp.inc` inside `clubb_kokkos`, before the generated
   metadata type. These are the three actual Native public scalar parameter
   types, each with eight ordered species fields/defaults. Root can instead
   teach the existing generator owner to emit those exact types in that order.
3. Keep the exact generated `hm_metadata_type` and its converters.
4. `default_tables.hpp.inc`, then `procedure_declarations.hpp.inc`, inside
   the same namespace. The latter includes thin generic View entry adapters.
5. Promote `corr_varnce_module.cpp`, add it to the core target, and replace
   matching lower call sites only after the reviewed dependencies are available.

The source currently expects those declarations to be merged; it is not a
standalone replacement header. No shared generator, live source or CMake file
was edited. `f2k_dependencies.hpp` declares one missing public lower contract:
`matrix_operations::mirror_lower_triangular_matrix(nvars,matrix)`, with required
integer `nvars` In and core-real rank-2 `matrix(nvars,nvars)` InOut. The two host
calls remain explicit, cloud then below. Generate that bridge through its
existing owner or port the canonical lower routine before promotion, then remove
the draft-only declaration include. Existing host-only
`f2k::clubb_at_least_debug_level_api(int)` is reused. No new selected constant,
metadata accessor, cache or device bridge is needed.

Metadata indices remain Native one-based, with existing negative/zero sentinel
semantics. Initialization resets the existing object, matching `INTENT(OUT)`
deallocation/default restoration, before saving all eight indices and the Ncn
variance ratio. It preserves the source's species-specific variance formulas,
name/flag/tolerance order and sequential PDF counting. Five local host mirrors
fill the generated View components; names are exactly ten characters including
Fortran blank padding. Valid Native input indices are unique and cover the
`hydromet_dim` enabled slots. Native leaves unassigned metadata slots undefined
for inconsistent/gapped inputs; Kokkos allocation may initialize those bytes.
Such invalid-input byte values are not a newly defined scientific API.

The two default correlation tables are immutable, Native-named `constexpr`
`std::array` values, as approved by root after checking every repository read
and write site. All 288 values retain Fortran-flat RESHAPE order; `_c` means
`core_rknd`, and all table/type literals retain their actual core-kind values.
The current standalone/JAX repository only reads these defaults. External
mutation of the former public Fortran arrays is not a supported new C++ API.
There is no mutable/thread-local C++ state or globally allocated Kokkos View.
Read-only unmanaged HostSpace matrix aliases exist only during the default call.

Correlation setup retains both-present optional override semantics, lower-to-upper
mirroring and the exact warning checks/messages. Small matrices stage locally
on the host and return into caller-owned output allocations. The symmetry check
preserves Native column/row order, asymmetry-before-diagonal precedence, fatal
assignment to all error columns, distinct header-array/global-header choices,
full first-index-fast matrix payload and first-failure return. It does not add
a finite-value check or change the Native tolerance. Diagnostic values use a
local stream at `max_digits10`; exact compiler-specific list-directed decimal
spacing is not claimed.

Matrix printing retains the Native `(m,n)` orientation, F5.2 fields stored in
LEN10, six-character insertion stride and LEN72 lines, including width-overflow
stars and trailing blanks. Its fixed finite-value formatting is implemented;
exact decimal rounding/non-finite spelling and the Native out-of-range substring
behavior beyond twelve displayed variables remain unverified. Source tables
have twelve variables; no wider-print behavior is claimed.

All substantive descriptions, equations, species comments and diagnostics remain
in Native order. Comments correct copied initialization overview index/type names
and the assertion description that claimed STOP, although Native sets fatal and
returns. The full Lebo et al. (2015) citation remains with the parameter types.
No JAX-only compatibility metadata, lookup helper, boolean-only error interface,
new scientific helper or speculative FMA/signed-zero adaptation is introduced.

`SOURCE_OUTLINE.json`, `SOURCE_AUDIT.json` and `CLEANUP_AUDIT.json` separate source
and cleanup review from runtime evidence. Optional-parameter combinations,
zero/all hydrometeors, reordered species, warning/fatal branches, matrix display,
host/strided callers, CPU/GPU execution and model equality are untested. Parent
integration owns the maintained CPU zero-tolerance model gate and GPU milestones.
