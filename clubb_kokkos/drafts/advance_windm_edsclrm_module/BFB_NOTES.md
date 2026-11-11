# Arithmetic and verification notes

This isolated draft has no compilation, runtime or bit-for-bit evidence. Root
must pass CPU Fortran versus CPU Kokkos with zero tolerance before integration
advances. GPU comparisons are occasional milestone evidence, with any differing
results documented separately. Source audits cannot establish either result.

No compiler-specific FMA, forced-rounding or reassociation correction was added.
The draft retains Native expression grouping, explicit Coriolis term order,
Crank-Nicholson halves, nudging order and direction-dependent RHS expression
order. Integer squares use explicit repeated-operand multiplication; core-kind
literals use `real_t`. The wind-speed maximum retains the positive Native `eps`
threshold. Its ordinary conditional lowering chooses the second operand on ties;
no extra signed-zero correction is introduced.

The Native direction-dependent RHS comments explicitly describe matching both
grid orientations at `-O0`. The draft preserves those equations. Native optional
band/row reversal is retained only at the existing solver storage boundary,
with its original side effects and fatal-return placement. No floating-point
reduction is introduced by the per-column status validators.

The retained Native/JAX definedness and scalar-flux distinctions are listed in
[INTEGRATION.md](INTEGRATION.md). No value is invented for unwritten Native
scratch. These source choices need model-level verification when reached; they
are not grounds for assuming BFB or changing the Native physics.

[source_audit.json](source_audit.json) records current source/header hashes and
contract/device/diagnostic review.
[formula_audit.json](formula_audit.json) separately records operation trees,
affine array addresses and ordered predicates. Its fixed three-band LHS slice
is compared as the equivalent explicit Native `j=1..3` loop; the actual solver
packing/reversal is reviewed manually. No numerical test suite was run for this
draft, as required by the independent-draft task.
