# Mean-advection operator draft

This is the complete supported source draft of
[`mean_adv.F90`](../../../src/CLUBB_core/mean_adv.F90): both public routines,
thermodynamic centered/upwind choices and momentum centered differencing.
The full 519-line Native file and 324-line JAX counterpart were read before
implementation. No source option is omitted. No build or numerical test was run.

Read [INTEGRATION.md](INTEGRATION.md) before promotion. Root owns the live-call
hookup, obsolete lower bridge removal and regular CPU Fortran versus CPU Kokkos
bit-for-bit gate.

Separate source review checks the ten- and seven-argument signatures/public
declarations, 32 assignment operation trees, four ordered predicates, ten local
diagonal/weight constants and all nine device bodies. Native numeric boundary
order, weight lower bounds, complete output ownership and comments were also
reviewed. There are no lower science calls, stats calls or explicit error paths.

Artifacts: [SOURCE_OUTLINE.json](SOURCE_OUTLINE.json),
[source_audit.json](source_audit.json), [formula_audit.json](formula_audit.json)
and [cleanup_audit.json](cleanup_audit.json).
[prepare.py](prepare.py) prepares this isolated draft;
[audit_source.py](audit_source.py) and [audit_formulas.py](audit_formulas.py)
compare source only. Arithmetic choices and pending numerical evidence are
recorded in [BFB_NOTES.md](BFB_NOTES.md).
