# Diffusion operator draft

This is the complete supported source draft of
[`diffusion.F90`](../../../src/CLUBB_core/diffusion.F90): both public routines,
both diffusion-option branches and all boundary rows. The complete 1037-line
Native file and 747-line JAX counterpart were read before implementation.
No source option is omitted. Compilation and numerical verification are pending.

Read [INTEGRATION.md](INTEGRATION.md) before promotion. Root owns the live-call
hookup, missing selected flag, obsolete bridge removal and regular CPU Fortran
versus CPU Kokkos bit-for-bit gate.

Separate source review checks both ten-argument signatures/public declarations,
59 assignment operation trees, seven ordered predicates, both five-argument grid
derivative calls, six private diagonal constants and all 23 device bodies.
Native numerical diagrams/comments, band ownership, unconditional momentum
scratch work and generic array entry overloads were reviewed separately.
No build, model execution or numerical test was run.

Artifacts: [SOURCE_OUTLINE.json](SOURCE_OUTLINE.json),
[source_audit.json](source_audit.json), [formula_audit.json](formula_audit.json)
and [cleanup_audit.json](cleanup_audit.json).
[prepare.py](prepare.py) prepares the isolated source;
[audit_source.py](audit_source.py) and [audit_formulas.py](audit_formulas.py)
perform source comparisons only. [BFB_NOTES.md](BFB_NOTES.md) states the pending
numerical evidence and arithmetic choices.
