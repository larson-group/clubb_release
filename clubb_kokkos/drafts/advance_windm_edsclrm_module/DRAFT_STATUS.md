# Wind and eddy-scalar draft

This is the complete supported source draft of
[`advance_windm_edsclrm_module.F90`](../../../src/CLUBB_core/advance_windm_edsclrm_module.F90).
It has not been compiled, integrated or numerically tested. Root owns promotion
and the required CPU Fortran versus CPU Kokkos bit-for-bit comparison.

Read [INTEGRATION.md](INTEGRATION.md) for missing public lower definitions and
constant owners. The source and header retain all six Native routines in order,
with all 52 required arguments of `advance_windm_edsclrm`. The five private
helpers remain in the unnamed namespace. No whole routine delegates to Fortran.
Only the UV sponge block is excluded, matching JAX scope and the existing
upstream C++ driver rejection.

The separate source review checks routine/signature order, 61 ordered calls,
427 non-solver argument mappings, both complete solver mappings, 115 explicit
and seven whole-array assignment operation trees, 58 ordered IF predicates and
39 diagnostic write payloads with all 72 value tokens. All 78 device bodies were checked
for owning-UDT capture, host access, backend calls and nested launches. Source
review and 100-column formatting do not establish compilation or numerical
agreement. The solver storage mapping and its error/side-effect order were
also reviewed manually against the complete Native body.

Review artifacts are [SOURCE_OUTLINE.json](SOURCE_OUTLINE.json),
[source_audit.json](source_audit.json) and [formula_audit.json](formula_audit.json). The separate cleanup record is
[cleanup_audit.json](cleanup_audit.json).
[prepare.py](prepare.py) prepares this isolated draft;
[dependency_contracts.py](dependency_contracts.py) records actual Native public
lower contracts. [audit_source.py](audit_source.py) and
[audit_formulas.py](audit_formulas.py) perform source comparisons only. They do
not run model physics. Numerical status and arithmetic choices are recorded in
[BFB_NOTES.md](BFB_NOTES.md).
