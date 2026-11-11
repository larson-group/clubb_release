# Draft arithmetic and BFB status

No explicit FMA, forced arithmetic association, alternate library or new floor
was introduced to repair a numerical difference. There has been no runtime
comparison of this draft, so its arithmetic/BFB status is pending integration.

Natural translations retain complete integer-power bases with parenthesized
square/cube products. Real powers, including `F_w**three_halves` and the unused
exponential blend's `abs(Skw)**lambda`, remain ordinary `Kokkos::pow` calls.
No speculative `sqrt`, reciprocal, `cbrt`, `log1p` or `expm1` rewrite is used.
Every surrounding sum, product and division group is separately source-audited.

The deliberate intrinsic mapping follows the root's earlier grid evidence:
Native MIN/MAX selects its second argument on an equality tie. Ordered C++
ternaries retain that choice at the responder variance zero floors, the passive
scalar correlation maximum, the maximum over all five responder correlations,
the correlation-flux cap, and the nonzero-skewness minimum `F_w` floor.
The five-argument MAX is folded in source operand order using its existing
`max_corr_w_x_sqd` variable. A local scalar candidate evaluates the passive-
scalar MAX's quotient only once, as Native evaluates an intrinsic argument.
These are explicit intrinsic/language boundaries, not newly demonstrated
compiler-rounding fixes. No new numerical experiment was run for them.

The preserved concrete MAX failure is
`/tmp/clubb_core_lower_validation/grid_operators_before_zero_tie.log`
(`+0` versus Native `-0`); its source review is
`/tmp/clubb_core_lower_validation/grid_source_audit.json`. This draft carries
forward that established comparison convention without independently
validating every new intrinsic call.

The integer `ANY` reduction contains no floating-point accumulation. Independent
array assignments share kernels; coupled arithmetic, source column order and
the driver fatal boundary retain Native order. Native implicit-coefficient
metadata and caller allocations are preserved.

If integration reveals a difference, root should record the exact expression,
first differing values, Native compiler/assembly evidence and maintained model
comparison result in `clubb_kokkos/BFB_ADAPTATIONS.md` before treating a special
arithmetic adjustment as verified.
