# Integration contracts

Read [DRAFT_STATUS.md](DRAFT_STATUS.md) first. This directory is a source draft,
with no compile or numerical evidence. Promote reviewed pieces into their
native-named owners only after the active CPU comparison is green. Root owns
production source, CMake changes, dependency deduplication and bridge removal.

The public driver in `new_pdf_main.hpp/.cpp` retains all 29 Native arguments,
including `nz`, `ngrdcol`, all required outputs, clipped `Skrt`/`Skthl` and the
caller-owned `implicit_coefs_terms`. Real fields use `Real2D(column,level)` and
LayoutLeft. The nine coefficient/explicit-term component Views are captured
individually; the owning derived type is never captured in a device closure.

`new_pdf.hpp` contains all nine Native lower procedures, including private
`sort_roots`. Public profile overloads retain Native argument order and required
`nz`. Independent-cell overloads omit only `nz`, use scalar reference outputs,
and are device callable. The limits routine currently has only a View/profile
API because its complex root dependency has not been promoted. Private main
helpers remain in an unnamed namespace and in Native order. C++'s reserved
`lambda` keyword is spelled `lambda_`; no argument is optional.

## Root dependency

The only unported lower science dependency in this draft is public
`calc_roots::cubic_solve`, called twice in `calc_limits_F_x_responder`:

```text
cubic_solve(nz, a_coef, b_coef, c_coef, d_coef) -> complex(nz,3)
```

Every argument is required: `nz` is integer IN, each coefficient is a
real(core_rknd) rank-one IN profile. Root order matters. Two limit branches
select the original first root's real part rather than the sorted largest root.

`f2k_dependencies.hpp` declares the current pending boundary only. Its four
inputs are owned LayoutLeft HostSpace real profiles and its return is an owned
LayoutLeft HostSpace View of `std::complex<real_t>` with shape `(nz,3)`. It has no
implementation, global storage, parser extension or fallback. Coefficients are
computed on device, packed through the existing `to_host_space` owner, and only
real parts of the returned roots are copied back. Each solve is a whole-profile
call; no parcel-by-parcel native crossing is introduced. Root may skip creating
this transient adapter by promoting the C++ root API directly.

Root's separate [`calc_roots` draft](../calc_roots/calc_roots.hpp) now publishes:

- `cubic_solve(nz, a_coef, b_coef, c_coef, d_coef)` for rank-one Views returns
  owned `Kokkos::View<Kokkos::complex<real_t>**,LayoutLeft,A::memory_space>`.
- Its device scalar overload returns `Kokkos::Array<complex_t,3>` for scalar
  coefficients in the same order.

At promotion, replace the two explicit calls with that public API, remove the
HostSpace coefficient/root staging, retain `.real()` extraction and the Native
sort/selection order, and remove this draft-only dependency declaration.
Do not edit the root-owned draft here or assume it already has runtime evidence.

## Constants, state and error behavior

Existing generated constants and parameter indices suffice. Tunable indices are
already zero-based in `parameter_indices.hpp`; do not subtract one again. The
main host orchestration mirrors `clubb_params` once because Native passes four
individual tunable values as scalar inputs to profile helpers. Other fields
remain on the execution space apart from the explicit root boundary.

All three model flags are required. The current selected generator exposes
`l_explicit_turbulent_adv_xpyp`; it must also expose Native
`l_explicit_turbulent_adv_wp3=false` and `l_explicit_turbulent_adv_wpxp=false`
through the existing generator owner, as already needed by the hybrid draft.
Both explicit and implicit branches are present in this draft.

Invalid `F_x`/`Skx` states retain Native fatal behavior. The profile mixture API
uses an integer validity reduction, prints both original lines once on the host
and throws at the fatal boundary before returning any mixture field. This check
precedes setter component calculations. Device scalar calls use `Kokkos::abort`
with both Native messages. JAX's nonfatal callback and symmetric substitute for
an invalid skewed state are not copied. C++ exceptions are a runner boundary
adaptation of Native `ERROR STOP`.

There are no statistics-backend calls in these two Native modules. `F_w`, `F_rt`,
`F_thl` and their limits retain their declared diagnostic outputs/scratch despite
Native noting that their output functionality is deprecated. No new stats UDT
or direct netCDF dependency is needed.

## Ownership and cleanup review

All scratch is call local. Profile APIs stage through the shared
`array_types.hpp` owner, restore every Out/InOut allocation, and return function
results in the first array input's memory space. They accept strided column
subviews and host Views without exposing device storage to a host result caller.
The parent retains Native column/call order; independent array assignments in the
same phase are combined where straightforward. Profile kernels are provisional
execution mappings, with no measured GPU performance claim.

The separate audit compares Native/target signatures and order, all 122
scientific assignment trees and predicates, every parent/helper call list,
complex-root coefficient mappings, original versus sorted-root selection and
all nine packed derived-type fields. View indexing/return adapters and device
scalar forms are necessary language/execution boundaries; no general helper
framework, silent fallback or scientific option removal was introduced.
