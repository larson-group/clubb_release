# Polynomial root draft

[calc_roots.hpp](calc_roots.hpp) implements the complete three-routine
[Native module](../../../src/CLUBB_core/calc_roots.F90). Source review and
formatting are recorded in [source_audit.json](source_audit.json); compilation,
runtime, BFB and performance evidence are pending. No production or shared
build/generator files were changed.

Promote the header into the core source directory and include it at the two
pending whole-profile root boundaries in the ADG/Luhar and new-PDF drafts.
`cubic_solve(nz,a,b,c,d)` returns owned `(nz,3)` roots;
`quadratic_solve(nz,a,b,c)` returns `(nz,2)`. Coefficients are rank1 Views;
results retain the first coefficient View's memory space and LayoutLeft.
All coefficient Views must be accessible from that View's execution space.

Roots use `Kokkos::complex<real_t>` on host and device. Replace the old draft
host-root `std::complex` aliases with that value type and use `.real()`,
`.imag()` and `Kokkos::abs` at their scalar inspection sites. Remove the
declaration-only temporary complex f2k interfaces when their final callsites
are replaced. No split-complex ABI conversion or shared parser extension is
needed for the direct C++ calls.

Device scalar overloads take four cubic or three quadratic real coefficients
and return `Kokkos::Array<complex_t,3>` or `<complex_t,2>`. They expose the
Native independent per-level equations for use inside another kernel without
launching a nested profile operation. The original profile API definitions
retain Native source order, followed by `cube_root` and the scalar execution
overloads. This is the only new execution adapter; no alternate root algorithm
or target-only polynomial guard was introduced.

Current Native already has the S/T sign correction needed for negative real
Cardano arguments. The draft retains `SIGN` through `copysign`, complex square
roots and a complex one-third exponent. It does not import JAX's safe-root or
gradient helper. Nonzero leading coefficients are a Native precondition.
Complex math uses ordinary Kokkos functions; no BFB-specific arithmetic change
has been added. Root ordering and complex-library rounding remain untested
until promotion and the relevant model/edge checks.
