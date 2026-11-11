# TSDADG PDF draft

[`new_tsdadg_pdf.hpp`](new_tsdadg_pdf.hpp) ports the complete supported
[`new_tsdadg_pdf.F90`](../../../src/CLUBB_core/new_tsdadg_pdf.F90) module. All four
routines retain native order, argument names, required status and field order:

| Routine | Visibility | Arguments |
| --- | --- | ---: |
| tsdadg_pdf_driver | Public profile API | 25 |
| calc_setter_parameters | Public device scalar API | 13 |
| calc_L_x_Skx_fnc | Public device scalar API | 6 |
| calc_respnder_parameters | Native private device scalar helper | 12 |

The native `respnder` spelling is preserved. The private helper uses the same
named detail-namespace pattern as the existing saturation and new-PDF headers;
its definition remains after the public routines. Output scalars are references,
and the public driver accepts the existing rank-one column subviews. There are
no new optional arguments, model flags or lower physics dependencies. Native
does not call a polynomial-root routine here, so no root bridge or alternate
root algorithm is introduced.

These files are isolated from the build. Promote the reviewed header into the
normal core directory, include it in `pdf_closure_module.cpp`, replace its
per-column `f2k::tsdadg_pdf_driver` call with the same direct C++ call, and remove
the obsolete whole-driver bridge through the current generator owner. Existing
`array_types.hpp`, `constants_clubb.hpp` and Kokkos math supply all dependencies;
the header-only implementation needs no new CMake source entry.

## Source behavior and execution

The driver retains the native profile dimension `nz`. Caller fields keep their
column-first LayoutLeft owner; generic profile Views accept strided column
aliases. The input View selects the execution and memory spaces. Inputs and
outputs must be conformable and accessible from that execution space.

The setter is chosen using native branch order: w wins ties with rt or thl;
otherwise rt wins a tie with thl; otherwise thl sets the mixture. Each branch
calls its setter first, then both responders in source order. The three
mean-spread calls, both covariance-sign branches, all four setter epsilon
branches, signed denominator floors and both component-variance formulas remain
explicit. Direct square roots and the original repeated formula associations
are retained; no JAX gradient-safe square root, eager evaluation of all branches,
new floor, explicit FMA or forced rounding is imported.

The native per-level signs, spread parameters and variance coefficients are
consumed only at their own level. They become distinct work-item scalars in one
profile kernel. The six immutable `small_l` values remain call-local constants
at `0.75`/`0.5`. Public scalar helpers allocate nothing, launch no nested kernels
and call no Fortran. Means retain their variable units; variances have squared
units, and normalized quantities are dimensionless.

All six negative-variance checks retain native order: w1, w2, rt1, rt2, thl1,
thl2. Each check saves its original negative variance, zeros that variance and
its coefficient, and increments an integer warning count. NaN and signed-zero
comparisons retain the native `< 0` condition. When the count is nonzero, a local
`(level,6)` LayoutLeft payload is copied to the host and the exact warning/data
labels print in native level/field order with the pre-clipping value. This host
I/O boundary replaces native printing inside the level loop; no warning is
silently dropped. There is no native statistics call or error-status object in
this module and no new global state.

Both full mathematical-method comment sections remain in their original source
positions. Two copied procedure-description errors are corrected: the setter
returns variances, and the responder consumes an existing mixture rather than
returning a newly set mixture. Their formulas and signatures are unchanged.

## Source review and pending numerical checks

[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json) records the complete 931-line native
read and actual JAX scope before implementation. A separate final
[source audit](SOURCE_AUDIT.json) and [expression audit](expression_audit.json)
cover four signatures, 53 ordered assignment trees plus six immutable constants,
14 native IF predicates, 12 helper calls with 129 arguments, all 12 warning/data
labels, six original-value payload mappings and both complete derivations.
The final cleanup preserves source argument groups, comments and readable
arithmetic associations with a 100-column limit.

No build, unit test, model run, BFB or performance evidence is claimed. All
setter/tie, near-epsilon, covariance-sign, variance-warning and memory-space
paths remain numerically untested. Default ADG1 model runs do not exercise this
TSDADG branch. After promotion, select the existing TSDADG PDF option in the
maintained CPU Fortran/CPU Kokkos model comparison and retain zero tolerances.
GPU comparisons remain occasional milestones. Record any later deliberate
arithmetic adaptation in the existing [BFB log](../../BFB_ADAPTATIONS.md).
