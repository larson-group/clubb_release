# Saturation draft status

The full native saturation module is translated in `saturation.hpp`, including
all ten native implementations and the four public generic/scalar APIs. All
scalar APIs are device callable. Liquid and ice arrays use column-first
`Kokkos::View` fields, execute on the default execution space and return owned
storage in the pressure field's memory space.

This is an isolated source draft. It has not been compiled, integrated or run.
No numerical, BFB or performance result is established. Production, CMake,
generators, tests and the root worklog were not changed by this task.

Completed source review:

- Read all 1185 native source lines and the actual JAX saturation implementation.
- Compared all ten native signatures and relative routine order; the sole native
  optional argument remains optional and 1-based.
- Compared 18 saturation/adjustment expression trees and nine mixing-ratio
  expression trees, including nested GFDL
  powers and the factored/Horner polynomials. Compared all 156 lookup literal
  tokens exactly and checked the shared scalar/array liquid vapor-pressure body.
- Re-read selectors, minimum temperatures, pressure guards, specific-humidity
  formulas, adjustment iteration/error semantics and scientific comments.
- Applied the repository C++ formatter and checked delimiter balance, whitespace
  and 100-column lines. These are source checks, not compilation evidence.

`source_audit.json` records source/header hashes and the exact source-only
checks. `audit_source.py` recreates that record without executing physics.
`prepare.py` recreates the first draft from canonical source and must be followed
by formatting and re-audit; do not rerun it over numerical fixes made later.
`INTEGRATION.md` lists owner dependencies and public contracts. `BFB_NOTES.md`
separates ordinary translations, source differences and pending numerical work.

All selectors, the currently disabled specific-humidity branches, both grid
directions, partial/empty ranges, lookup rounding edges and adjustment failure
paths remain runtime untested. Root owns incremental integration and the
maintained whole-model Kokkos/Fortran comparison.
