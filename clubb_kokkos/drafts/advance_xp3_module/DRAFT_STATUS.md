# Third-moment draft status

The complete supported `advance_xp3_module.F90` is translated in the native-named
C++ source/header. Both public parents and all three private routines are present
in native order. No native scientific option is excluded in this module.

The separate source review covers all five signatures (38, 28, 18, 7 and 4
arguments), all 69 calls including grid/scalar functions, all 425 call arguments,
nine central expression trees and nine kernel bodies. ADG1/non-ADG1 scalar
minima, the predictive/LMM and steady-state branches, original stats names/gates,
last-storage-level boundary and full scientific comments were reviewed separately.
The repository formatter, delimiter balance, whitespace and 100-column checks
pass. `source_audit.json` records current hashes and exact source-only evidence.

No compilation, numerical tests or timing runs were performed. This is not an
integrated or BFB-validated port. All branches remain runtime untested here,
including the source's compile-time-disabled predictive branch. Production
source, CMake, generators, tests and the root worklog were not edited.

Existing live grid, Skx, stats, array and constant owners cover every lower
requirement. There are no missing bridges, speculative declarations, new
constants or new state owners. Root owns promotion and regular CPU Fortran/CPU
Kokkos zero-tolerance model checks, with occasional GPU checkpoints.

`SOURCE_OUTLINE.md` and `native_inventory.json` were written before implementation.
`INTEGRATION.md` lists actual entry points and necessary execution adaptations.
`BFB_NOTES.md` records ordinary arithmetic choices and pending numerical evidence.
`audit_source.py` audits source structure without executing physics. `prepare.py`
recreates this initial draft; it must not overwrite later integration fixes.
