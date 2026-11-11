This draft implements the three Native model-flag procedures, with a same-named
stream overload for the actual C++ driver's output sink. The complete 1,040-line
Native module and 322-line JAX behavior mirror were read before the outline.
JAX's initializer is implicit in `ConfigFlags` construction; the C++ port retains
the original procedure and both Native fatal selector checks.

Native has no mutable module flag instance: it declares immutable constants and
the configuration type. Existing `ClubbDriver` scalar flags and configuration
object already own mutable case state. The default routine writes all 66 scalar
outputs in Native order; the initializer checks grid-remap and mixing-length
selectors before assigning all 66 fields. Native invalid IntentOut data are
undefined after a fatal stop; the C++ port does not fabricate initialization for
that invalid path. No additional state class, mutable global or duplicate UDT
is introduced.

Merge the declaration/include fragments beside the authoritative generated type
and compile `model_flags.cpp` in the permanent core library. Extend the existing
selected-constant producer with `ppm_remap` and `no_grid_adaptation`; their source
values are 2 and 0. Preserve the full generated constants-and-type region using
the existing shared marker design. Its begin marker must precede the constants
loop, so future regeneration updates constants instead of preserving stale prefix
bytes. Permanent APIs/includes remain outside that region.

The canonical integer-`iunit` printing contract remains an explicit public Native
I/O service, via the staged existing-manifest entry. A `std::ostream& iunit`
overload implements the source row/value sequence for the existing C++ driver's
`std::ofstream`/stdout ownership. It uses call-local text storage and needs no
logical-unit registry. Native prints 64 rows: it omits `saturation_formula` and
`fill_holes_type`, and its `l_mono_flux_lim_um`/`vm` labels use the opposite wind
field. These actual source writes are retained; no Native correction is made.
The stream adapter emits decimal integers and T/F logical values without changing
caller stream formatting; compiler-specific list-directed padding is unverified.

At hookup, replace the driver's whole default/initializer bridges directly.
Its current manual 66-row dump differs from Native; replace that duplicated dump
with the stream overload at the existing stdout/case-log boundary. The existing
Native configurable-flag namelist reader remains a separate I/O service.
Retire the two obsolete whole-procedure bridge entries after their final callers
move, retaining only the genuine integer-unit printer where needed.

`SOURCE_AUDIT.json` checks signatures, all 132 assignments, two fatal checks and
64 ordered print mappings. `CLEANUP_AUDIT.json` records a separate source-group,
comment and whitespace pass with unchanged executable tokens. No deliberate BFB
arithmetic adaptation was added. Independent review, compilation, all flag/error
branches and text-output runtime evidence remain pending group integration.
