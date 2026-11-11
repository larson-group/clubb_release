# Polynomial root source outline

Read the complete316-line Native `src/CLUBB_core/calc_roots.F90` and its JAX
counterpart before implementation. Three public routines, no private Native
helpers, stats, diagnostics, saved state or optional arguments:

1. `cubic_solve(nz,a_coef,b_coef,c_coef,d_coef)` returns complex `(nz,3)`.
   Preserve Q, R, discriminant, complex square root, signed S/T complex cube
   powers and the three roots in exact source order. The current Native source
   already contains the sign correction for negative real Cardano arguments;
   do not copy JAX's separate safe-root/gradient helper or reorder roots.
2. `quadratic_solve(nz,a_coef,b_coef,c_coef)` returns complex `(nz,2)`.
   Preserve discriminant, complex square root and plus/minus root order.
3. `cube_root(x)` is scalar real, using `pow(x,one_third)` for x>=0 and
   `-pow(abs(x),one_third)` otherwise. Preserve negative-zero comparison and
   the real-power operation; do not replace it with a different cbrt algorithm.

Profile coefficients are rank1; the first return dimension is vertical level,
as Native has no column dimension. Outputs retain the input memory space and
LayoutLeft. Native requires nonzero leading coefficients and conformable
arrays; introduce no new coefficient clipping, degeneracy guard or tolerance.

Required language adaptation: expose device scalar cubic/quadratic overloads
so callers inside Kokkos kernels can evaluate one polynomial without launching
a nested profile kernel. Native-signature profile definitions remain first in
source order; scalar overloads at the end own the unchanged per-level equations.
Forward declarations resolve C++ lookup. Results use Kokkos::complex<real_t>
in both host and device memory; older draft std::complex host adapters need an
explicit type/callsite update at promotion, not a new Fortran ABI shim.

All numerical and build evidence is pending. This is an isolated draft; no
production, shared CMake/generator or bridge changes are authorized by this
outline alone. No speculative FMA or special math workaround is planned.
