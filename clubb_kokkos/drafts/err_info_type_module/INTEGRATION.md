This draft implements the five canonical error-context lifecycle routines in
Native source order. The two Native private-context converter accessors remain
transient bridge support and are excluded from the permanent C++ API. The full
415-line Native file and actual JAX `err_info.py` were read before the outline.
JAX's MPI default, on-demand headers and omitted coordinate normalization do not
replace the Native contract.

Merge `procedure_includes.hpp.inc` at the top level of the existing generated
`err_info_type_module.hpp`, then merge `procedure_declarations.hpp.inc` inside
`clubb_kokkos` after its authoritative mirror type. Replace the old generated
inline `cleanup_err_info_api`; do not keep both definitions. Compile this folder's
`err_info_type_module.cpp` as the permanent procedure owner. Root must preserve
these regions in the existing type generator rather than duplicate the schema.

Initialization and optional state updates run in C++. Default MPI/chunk values
are both one, coordinates use `unused_var`, and allocation/reset/update order
matches Native. Required initialization accepts only unallocated component
storage: an already allocated object throws rather than silently reallocating.
Native partial-allocation failures may modify earlier components before stopping;
that invalid-input failure progression is not reproduced. Cleanup releases Views
in Native order but has no Fortran DEALLOCATE STAT status or corresponding warning.
Scalar context and the global header string remain unchanged on cleanup.

The setter retains Native's actual repeated `err_header` allocation test, rather
than adding checks for `host_lon` or `err_code`. The optional updates remain
MPI, chunk, latitude, longitude. The exact `>90`/subtract-90 and `>180`/subtract-180
branches are preserved, including strict ties and unmodified other values.

Exact fixed-width header I/O temporarily remains in Fortran. Merge the explicit
alias from `f2k/procedure_manifest.py.inc` into the existing generator manifest:
`format_err_info_headers` selects PUBLIC `set_err_info_values_api` with an empty
optional selection. Add its alias-to-public-name entry to existing `PUBLIC_NAMES`,
just like the stats presence overload. Both generated C++ and Fortran retain
`set_err_info_values_api`; the distinct C symbol still selects the empty optional
set. Private C++ `update_headers` calls its two-required-argument f2k overload.
No new alias mechanism is needed. The original setter's PRESENT guards then skip all state
updates and normalization, while its unconditional private `update_headers` call
formats the copied case context. C++ private `update_headers` invokes this service
using a call-local snapshot and returns only column/global header data. The
original column-header allocation is retained. No private Native formatter is
exposed, no new private accessor is added, and no mutable global state is created.
The existing derived-type converter retains all 200/214 character bytes, including
padding. Native assumes at least one column when formatting the global header.

Promote existing driver initialization calls by dropping `f2k::`. The current
latitude/longitude-only setter call needs two `std::nullopt` placeholders for the
Native optional integer positions before its two Views. Remove the old whole-init
and selected lat/lon setter manifest entries only after their last callers move;
keep the temporary formatter alias and the converter's existing private accessors.
They can disappear together when fixed-width formatting is ported.

`SOURCE_AUDIT.json` covers all signatures, defaults, allocation/cleanup order,
optional branches and both coordinate kernels. `CLEANUP_AUDIT.json` records the
separate source-group/comment/whitespace pass and unchanged executable tokens,
with unused-include removals listed separately. No arithmetic contraction or
other deliberate BFB adaptation was introduced. Independent review, compilation
and all lifecycle/error/format branches remain untested in this draft.
