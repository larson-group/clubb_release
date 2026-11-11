# Native metadata lookups

All five `index_mapping.F90` functions are implemented in Native order. Promote
`index_mapping.cpp` and `.hpp`, add the source to the core target and replace
matching explicit lower boundaries in the reviewed callers. The existing
`corr_varnce_module.hpp::hm_metadata_type` and its converters remain unchanged.
All radius constants and `zero` already exist in `constants_clubb.hpp`; no new
bridge, constant, metadata type, global cache or accessor is required.

These functions run on the host and read only scalar metadata. The current
fill-holes draft precomputes the two needed lookups before its level kernel;
replace those call sites directly, without moving an owning metadata object
into a GPU capture. Other callers must likewise keep this host metadata boundary
unless a later measured device requirement establishes a distinct safe API.

Indices remain Native one-based. Unmatched lookups return zero; comparisons
against absent metadata still follow the original ordered IF/ELSEIF chain.
There is no additional positive-index guard. The inverse functions deliberately
use different Native species orders. JAX's guarded zero-based tuples/helpers
are secondary scope evidence and were not imported. The maximum radius result
retains metres, the original species OR predicates and zero initialization.

Full scientific descriptions and species comments remain. The old module
comments naming `array_index`, `setup_pdf_indices` and `init_microphys` are
corrected to the actual `hm_metadata_type`/initialization owners inspected in
the Native source. No executable behavior is changed by those corrections.

`SOURCE_OUTLINE.json`, `SOURCE_AUDIT.json` and `CLEANUP_AUDIT.json` record the
source review. Compilation, all lookups/sentinel combinations, CPU/GPU caller
execution and whole-model equivalence are pending. No numerical test or BFB
claim is made, and no speculative arithmetic adaptation is present.
