This source draft has no compilation, runtime or numerical comparison evidence.
CPU Fortran versus CPU Kokkos BFB is the integration gate. GPU comparisons are
occasional milestone evidence under the accepted project policy.

No speculative FMA, forced association, floor, tolerance, reassociation or
alternative finite predicate was added. The only new real-valued calculations
are Native `rvm = rtm - rcm` and the exact left-associated
`(integral_after-integral_before)/dt + flux_top - flux_sfc - integral_forcing`
conservation residual. The settings reduction preserves the exact C1/C14
absolute-difference test and each `> one OR < one` comparison. Source audit
checks their expression trees independently.

Integer reductions replace logical ANY without changing boolean results. The
four existing live bodies are byte-identical to the starting snapshots.
Native all-column error propagation, diagnostic order, active and commented-out
fatal updates, and return placement are retained. No source-definedness repair
is silently applied. Inherited diagnostic quirks are explicit in
[`INTEGRATION.md`](INTEGRATION.md).

Diagnostic payload equality is a source claim. Host streams use C++ formatting
rather than Fortran list-directed number/boolean formatting; this is distinct
from saved numerical-field BFB. The host GPU-input fatal stop uses the existing
C++ exception mechanism at its Native compile-gated location.
