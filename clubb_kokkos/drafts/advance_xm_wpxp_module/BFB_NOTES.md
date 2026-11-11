# Mean/flux arithmetic and native definedness

No forced FMA, reassociation, alternate math approximation, new denominator
threshold or other speculative BFB workaround is present. This source draft has
not been compiled or compared numerically. The regular gate is CPU Fortran versus
CPU Kokkos at zero tolerance; GPU differences are recorded at occasional GPU
checks under the latest user instruction.

Ordinary translations preserve native formula trees: integer squares become
operand multiplication; exp/abs use Kokkos math; real literals cast to `real_t`;
ordered two-operand MIN/MAX use comparisons/ternaries with source operand order.
No cross-compiler NaN result is claimed for those intrinsics. Tunable-parameter
indices use the established zero-based owner; other native loop/grid indices
are converted only at View access or the existing bridge boundary.

Native over-implicit weights, first-pressure coefficients, field-specific
thresholds and argument choices remain literal. For example, thlm's clipping/
stats call still receives the exact lhs_pr1 field chosen in native source;
it was not harmonized by analogy with a neighboring call. The 328 independently
compared assignment trees and 1744 checked call arguments establish source
structure, not floating-point BFB results.

## Explicit source-definedness deviation

`xm_correction_wpxp_cl` writes `l_clipping_needed(i)=true` only for clipping
columns at native line 5917, then reads every column's mask at line 5945. There
is no false initialization. Its `xm_tndcy_wpxp_cl` output is also unwritten for
nonclipping columns before stats output at line 5958. The JAX implementation
instead derives a complete boolean column mask and sets unselected tendency
contributions to zero (advance_xm_wpxp_module.py:2075–2088).

This draft defines the complete mask false/true with an ordered per-column
scan, then defines nonclipping-column stats contributions as zero, matching
JAX. It retains the native early return when no column needs clipping and the
same physical masked formulas for columns that do. This is a documented
source-definedness correction, not a measured BFB workaround. There is no native
BFB claim for the previously undefined mask or tendency entries. Root owns any
native-source/backlog decision.

## Other inherited native edges retained

- In the compile-time-disabled explicit scalar TA path, native zeroing writes
  `lhs_ta_wpsclrp(...,sclr_dim)` inside the scalar loop. Its subsequent public
  RHS call receives rank-3 LHS storage rather than the scalar RHS. The draft
  retains that mapping with explicit native sequence association and leaves
  scalar RHS unwritten, while actual JAX assigns the RHS. Do not claim BFB for
  the source-unwritten RHS; changing that behavior requires a separate native
  decision. `l_explicit_turbulent_adv_wpxp` currently remains false.
- `xm_wpxp_solve.old_solution` is declared as `(nsup+nsub+1,ngrdcol,2*nzm-1)`
  at native lines 4330–4331, while its callers allocate/pack
  `(ngrdcol,2*nzm-1,nrhs)`. Only the optional descending BiCGStab branch indexes
  it directly. The draft keeps caller storage for the existing typed solver;
  that solver rejects BiCGStab in accordance with its established JAX scope.
  There is no claim for the inconsistent native indexed branch.
- Positive-definite adjustment stats temporaries are native-unwritten if no
  adjustment is needed, even though that source branch can sample them.
  The default `l_pos_def=false` is retained and no silent zero fill was added.
- Native diagnostics, including old prognostic fields, full matrix/RHS dumps
  and complete field dumps, remain in source order. Diagnostic formatting uses
  max_digits10 without altering the global output stream's precision.

Matrix axis permutation, whole-array reverse snapshots, device reductions,
subviews and host I/O are execution/layout adaptations listed in INTEGRATION.md.
They are not measured numerical workarounds. Root must log any later deliberate
arithmetic fix with expression, observed cause and before/after evidence in
`clubb_kokkos/BFB_ADAPTATIONS.md`.
