# Draft arithmetic and BFB notes

No FMA, forced arithmetic association, alternate math library or new numerical floor/denominator guard has been introduced to repair observed differences. There has been no numerical run of this full draft. Do not describe it as Native BFB yet.

Deliberate intrinsic mapping to carry into the permanent catalog at integration: native MIN/MAX calls select their second operand on equality, expressed with ordered C++ ternaries. This covers the responder alpha clamp, ADG1 mixture clamp, ADG2 sigma cap, Luhar m floors, backsolve alpha/coefficient/root floors, and largest-root coefficient/real-root maxima. The native tie behavior was established in the root's prior grid investigation; no new experiment was run for these expressions.

Natural translations: integer powers2/3 use multiplication of the complete source operand with source grouping; fractional `abs(Skx)**one_third` uses ordinary `Kokkos::pow`, preserving the source exponent rather than copying JAX's cbrt/gradient guard. Square roots remain ordinary Kokkos sqrt. Native HUGE maps to numeric_limits<real_t>::max(), preserving the finite native value rather than JAX's infinity. The real cubic-root selection preserves Native's unguarded linear division.

Language-boundary temporaries evaluate the fractional-power MAX operand and max_cubic_root function exactly once. They do not force a floating-point rounding/association choice. Host profile staging is required by the pending native complex-root boundary and is not evidence of CPU/GPU BFB or performance.

The sign initialization in close_Luhar_pdf is a separately documented source-definedness deviation, explicitly accepted by root and copied from actual JAX. It is not an arithmetic BFB workaround; Native normalized outputs on that nonvarying branch are undefined. Physical mean/variance formulas are unchanged.

When a runtime difference is actually found, record its exact routine/expression, first differing values, native compiler/assembly evidence, tested configuration/case and why the chosen change resolves it in the permanent clubb_kokkos/BFB_ADAPTATIONS.md. Keep these draft hypotheses distinct from demonstrated fixes.
