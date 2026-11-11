# Arithmetic and pending BFB evidence

No forced FMA, alternate root algorithm, reciprocal rewrite, clipped square
root or other speculative BFB workaround was introduced. The draft has no
compiled or numerical evidence. Source-tree agreement does not establish BFB.
Root must compare CPU Fortran with CPU Kokkos at zero tolerance before advancing
integration; GPU runs are occasional milestone evidence with differences logged.

Ordinary translation choices retained for review:

- Native integer squares/cubes become fully grouped multiplication trees.
  General powers, including the tunable exponent and fifth power, use
  `Kokkos::pow`; exponential, absolute-value and square-root calls use the
  corresponding Kokkos function. No contraction or library substitution is
  presumed to match merely because the build is unoptimized.
- All sums, products, quotients, setter variance formulas and implicit moment
  equations retain their Native operation grouping. The independent source
  audit checks 122 scientific assignment trees, including each expanded-power
  base and its surrounding groups.
- Ordered MIN/MAX use comparison/ternary expressions, choosing the second
  operand on a tie, consistent with the existing port's Native-intrinsic
  boundary. This is not new measured evidence for these routines or NaNs.
  Existing signed-zero motivation is recorded in the canonical
  [BFB adaptation log](../../BFB_ADAPTATIONS.md).
- Native signs use `>= zero` branches. Do not replace them with a sign-bit or
  `copysign` operation: negative zero follows the nonnegative branch here.

Source behavior preserved rather than replaced with JAX numerical adaptations:

- Invalid mixture/skewness states are fatal; the zero-spread symmetric limit
  applies only when Native reaches its valid zero-skewness branch.
- Responders retain both signed skewness caps and the factor 0.99, then derive
  spread limits from the two original cubic equations and sorting branches.
- The one-real-root branches use the real part of the original first root;
  the other branches use the explicit Native sorting network. Root order,
  equality comparisons, operand order and the final zero/one floors matter.
- Each implicit/semi-implicit coefficient retains its own zero-component or
  zero-variance branch. No JAX denominator substitute or clipped-gradient helper
  is added to an inactive or active branch.
- Setter interpolation and responder interpolation retain Native separate
  products; no repeated-term precomputation changes their operation tree.

Temporary root staging, scalar/array overloads, integer validity reduction and
host exception/device abort are documented execution/language boundaries in
[INTEGRATION.md](INTEGRATION.md), not demonstrated BFB improvements. Any later
mismatch-driven expression change needs actual failing/passing evidence and a
reviewable entry in `clubb_kokkos/BFB_ADAPTATIONS.md` owned by root.
