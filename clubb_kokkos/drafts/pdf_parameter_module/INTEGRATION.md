This draft implements the seven JAX-supported Native PDF lifecycle, packing
and diagnostic routines in source order. The full 899-line Native source and
326-line actual JAX `pdf_params.py` were read before the outline. The existing
47-field PDF and 27-field implicit-coefficient mirrors remain authoritative.

Merge `procedure_includes.hpp.inc` at top level and `procedure_declarations.hpp.inc`
inside `clubb_kokkos` after the existing generated types. The latter also places
Native `num_pdf_params=47` in its original owner. Compile the `.cpp` as the permanent
procedure owner. Root must preserve these regions through the existing generator;
this draft does not duplicate or alter either type schema.

Merge `clubb_api_module.hpp.inc` into the actual C++ API owner. This is the existing
Native three-argument `init_pdf_params_api` wrapper (Native lines 2437–2455), which
forwards directly to `init_pdf_params`. The live driver calls that original API
name twice. Keeping the forwarder avoids renaming its callers or inventing a
second initialization owner.

Initializers implement IntentOut by resetting only the caller-owned object,
allocating every component in Native order, setting metadata, then invoking the
corresponding zero routine. The 19 rank-2 implicit fields always exist after
initialization; eight rank-3 scalar fields are allocated/reset only when
`sclr_dim>0`. Reset kernels preserve allocations and per-cell field order, with
column-first LayoutLeft storage and individual View captures.

Packing/unpacking retains the exact 47-slot mapping of column one. The packet
shape is `(nz,num_pdf_params)` with `(level,field)` meaning, not a column axis.
Both optional bounds are Native one-based indices, used only when both are
present; one optional alone is ignored. The slice length must equal `nz`, as
required by Native's conformable array assignment; the C++ boundary explicitly
throws for a nonconformable length. It does not clamp bounds or add a safe fallback.
Unpack changes only the selected levels of column one, preserving other columns,
metadata and untouched levels, unlike the JAX zeroed reconstruction. Host/strided
packet adapters reuse the existing array owner once at the call boundary.

Printing stages the 47 fields once and emits every column/field/value in Native
order, retaining all labels. C++ diagnostic formatting keeps full real precision;
compiler-specific Fortran list-directed whitespace/number spellings are unverified.
No scientific calculation or error status is changed by this I/O adaptation.

`copy_single_pdf_params_to_multi` and `copy_multi_pdf_params_to_single` have no JAX
bodies and no repository call sites: explicit name-plus-parenthesis searches find
only their two owner definitions and two API definitions (Native owner lines
757/831, API lines 2465/2539). Their Native exports remain untouched; an external
host API port is not inferred from those exports. They are explicitly excluded
from this supported C++ group rather than replaced with fabricated behavior.

Replace existing driver/PDF-closure initialization and zero calls by dropping
`f2k::`, using the staged API forwarder where its original name is required.
Retire corresponding whole-procedure bridge entries only after their final live
callers move. `SOURCE_AUDIT.json` checks all 47 PDF field mappings plus the 19+8
implicit fields, optional behavior, printer order and original API wrapper.
`CLEANUP_AUDIT.json` records the separate source-group/comment/whitespace pass and
unchanged executable tokens. No deliberate BFB arithmetic adaptation was added.
Independent review, compilation and all lifecycle, optional-range and diagnostic
branches remain pending compatible-group CPU zero-tolerance validation.
