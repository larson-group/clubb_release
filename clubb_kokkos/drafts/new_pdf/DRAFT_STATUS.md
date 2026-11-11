# Original PDF draft

The complete supported bodies of `new_pdf_main.F90` and `new_pdf.F90` are present
in this isolated directory: the 29-argument driver, its three private helpers,
and all nine lower procedures. Native routine order, required argument order,
clipping, explicit/implicit branches and diagnostics are retained. Native remains
the authority; JAX was consulted only to confirm scope.

**This draft has not been compiled, integrated or numerically tested.** No BFB or
performance result is established. Production source, CMake, bridge generators,
tests and the root worklog were not changed by this task. Root owns promotion
and the maintained whole-model comparison.

The lower limits procedure currently declares two pending whole-profile complex
root calls. There is no implementation or fallback in this directory. Root has
prepared a separate `calc_roots` draft with the replacement public View API;
[INTEGRATION.md](INTEGRATION.md) records both contracts. This explicit dependency
must be resolved before compiling or integrating the full parent.

Completed source review:

- Read all 3086 Native source lines, outlined all 13 contracts and branches before
  implementation, and inspected the JAX scope/deviations.
- Independently checked all profile signatures, scalar adaptations, Native
  definition order, 122 scientific assignment trees, all branch predicates,
  parent/helper call arguments, root coefficient/selection equations, sort order
  and the nine packed derived-type fields.
- Re-read error behavior, all three advection flags, caller-owned outputs,
  View memory boundaries, scientific comments and citations.
- Applied the C++ formatter and checked 100-column lines and trailing whitespace.
  These checks inspect source; they do not execute physics or establish BFB.

[source_outline.json](source_outline.json) is the pre-implementation inventory.
[source_audit.json](source_audit.json) records reviewed hashes and source-only
checks; `python3 clubb_kokkos/drafts/new_pdf/source_audit.py` reproduces it.
`prepare.py` recreates the draft before formatting and audit. Do not rerun it
over later numerical changes without reconciling those changes first.

Root must establish **CPU Fortran versus CPU Kokkos BFB at zero tolerance before
advancing integration**. Use the maintained short multicolumn comparison during
integration; GPU runs provide occasional milestone evidence, with any differences
documented. All numerical branches and error paths remain runtime unverified.
