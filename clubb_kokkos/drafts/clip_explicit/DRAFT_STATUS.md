# Explicit clipping draft

This is the complete supported source draft of
[`clip_explicit.F90`](../../../src/CLUBB_core/clip_explicit.F90): all six public
routines, all 19 public integer constants and the genuine optional scalar
`threshold_hi`. No source option is omitted. It has not been compiled,
integrated or numerically tested.

Start with [INTEGRATION.md](INTEGRATION.md) before promotion. Root owns shared
registration, Native bridge removal and the regular CPU Fortran versus CPU
Kokkos bit-for-bit gate. GPU comparisons remain occasional milestone evidence.

The separate source review checks all six signatures and public declarations,
28 ordered calls with 117 argument mappings, the array-valued Heaviside call,
33 Native assignment operation trees, 26 automatically ordered predicates and
nine manually audited variance predicates. Both variance assignment paths are
checked. Four complete diagnostic payloads contain 38 value tokens, and all
32 device bodies are free of owning-UDT capture, host I/O/backend calls and
nested launches. These checks establish source structure, not compilation or
numerical agreement.

Review artifacts are [SOURCE_OUTLINE.json](SOURCE_OUTLINE.json),
[source_audit.json](source_audit.json), [formula_audit.json](formula_audit.json)
and the separate [cleanup_audit.json](cleanup_audit.json).
[prepare.py](prepare.py) prepares this isolated draft;
[dependency_contracts.py](dependency_contracts.py) inventories its one unresolved
public lower contract. [audit_source.py](audit_source.py) and
[audit_formulas.py](audit_formulas.py) compare source only and run no physics.
Arithmetic and remaining numerical evidence are described in
[BFB_NOTES.md](BFB_NOTES.md).
