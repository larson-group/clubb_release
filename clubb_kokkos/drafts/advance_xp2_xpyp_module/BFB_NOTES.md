# Arithmetic and pending BFB evidence

This draft has no compiled or numerical evidence. CPU Fortran versus CPU Kokkos
at zero tolerance is root's integration gate. GPU comparisons are occasional
milestones, with differences and causes documented. Neither the source audit
nor an unoptimized build guarantees BFB.

No explicit FMA, alternate special-function approximation, reciprocal rewrite,
reassociated pressure/advection formula or speculative precision workaround was
introduced. Integer squares use grouped multiplication; real powers and square
roots use ordinary Kokkos functions. Full Fortran parentheses and field-specific
products/quotients remain intact, including the rain-evaporation formulas where
JAX factors common expressions differently.

The existing ordered MIN/MAX convention selects the second operand on a tie.
Native `SIGN(one,value)` maps to `Kokkos::copysign(one,value)`, including negative
zero. These are source/language translations; there is no new measured BFB claim
for this module. Any mismatch-driven later change belongs with exact evidence
in root's [BFB adaptation log](../../BFB_ADAPTATIONS.md).

The following execution changes preserve source ordering constraints:

- Single-LHS equality uses integer disagreement reduction with the original
  subtraction/addition, division by two and multiplication by eps order.
- Matrix axis conversion and Native direction reversals copy values without
  arithmetic. Solver equations remain in the existing solver owner. The fatal
  return precedes final output reversal, as in Native.
- The precipitation-fraction recurrence remains sequential per column; it is
  not parallelized over dependent levels. Other cell operations preserve their
  statement and boundary order.
- Budget expression buffers evaluate the same source expressions immediately
  before the same Native calls. Explicit `l_count_sample=false` is retained for
  three clipping corrections rather than silently using the default.

Source-definedness is not a reason to claim numerical agreement. The dormant
microphysics limiter's stats temporary and source-unhandled PDF coefficient
selections are not repaired by inventing defaults. The Native sponge block is
excluded with JAX's supported scope and requires upstream rejection; its absence
is documented separately from BFB. All numerical branches, solver modes, grid
orientations, stats paths and error paths remain runtime untested in this draft.
