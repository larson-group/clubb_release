# Matrix operations draft

This isolated header implements the two routines supported by the actual JAX
path, `Cholesky_factor` and `mirror_lower_triangular_matrix`, plus the two Native
routines needed by complete Cholesky diagnostics:
`symm_covar_matrix_2_corr_matrix` and private `Symm_matrix_eigenvalues`.
Definitions retain Native order; the private helper uses a named detail namespace
and forward declaration. Root approved this four-routine scope.

Four unused, JAX-omitted public routines remain excluded:
`row_mult_lower_tri_matrix`, `set_lower_triangular_matrix`,
`get_lower_triangular_matrix` and `print_lower_triangular_matrix`. The latter
accepts arbitrary Native logical I/O units. No C++ unit/file registry or substitute
for that API was introduced; current diagnostics print their Native payloads
directly and do not need it.

Promotion requires four new public typed entries in the existing
`CLUBB_core/lapack_interfaces.hpp` owner, backed by the existing transient
`f2k/lapack_bridge.hpp`/Fortran raw-library shim. `lapack_dependencies.hpp`
records the declarations and should disappear once the actual owner supplies
those entries. All buffers are contiguous host storage with first index fast;
character selectors truncate to their first character in the shim, following
the already live expert-solver pattern. Dispatch to s/d libraries according to
`real_t`, matching Native generic selection:

| Public entry | Ordered argument contract |
| --- | --- |
| `lapack_potrf` | `uplo` character In, `n` integer In, `a(lda,n)` core real InOut, `lda` integer In, `info` integer Out |
| `lapack_poequ` | `n` In, `a(lda,n)` In, `lda` In, `s(n)` Out, `scond` Out, `amax` Out, `info` Out |
| `lapack_laqsy` | `uplo` In, `n` In, `a(lda,n)` InOut, `lda` In, `s(n)` In, `scond` In, `amax` In, `equed` character Out |
| `lapack_syev` | `jobz` In, `uplo` In, `n` In, `a(lda,n)` InOut, `lda` In, `w(n)` Out, `work(lwork)` workspace, `lwork` In, `info` Out |

The Native wrappers have old-style no-INTENT assumed-size declarations; these
contracts follow their complete bodies and the called library interfaces.
No raw library symbol/ABI logic belongs in this physics header. Existing host
`f2k::clubb_at_least_debug_level_api` remains explicit. No selected physical
constant, new metadata type or mutable cache is required.

Cholesky owns all host staging/scratch per call and restores both caller output
Views. It retains Native poequ/conditional-laqsy order, character flags, ten
total potrf attempts, info/error branches and complete ordered debug matrices,
eigenvalues and correlations. It restores the scaled/original matrix before
each modification, scans the diagonal in source order and preserves
`d_smallest*d_coef*real(iter,core_rknd)`. On the last positive failure it still
modifies and returns the matrix, as Native does; no extra retry or new fatal
check is added. Native illegal-argument ERROR STOP becomes a C++ fatal exception.
The eigen helper keeps fixed `lwork=180` and source convergence diagnostics.

The two independent kernels retain the caller memory/execution space. Mirroring
changes upper entries only; covariance conversion zeroes its output then fills
only the lower triangle. Inputs/outputs to the latter must be accessible from
the output execution space. Both can accept strided Views; Cholesky uses the
existing shared array owner for host packing/restoration. No device f2k or
global Kokkos allocation occurs.

Full scientific comments, references and diagnostic labels/payloads are retained.
Numeric diagnostics use local stream width/precision settings; exact
compiler-dependent G/E field spacing/rounding is unverified. No speculative
FMA, reassociation or signed-zero choice is present. Compilation, LAPACK status
branches, strided/host/device callers and all CPU/GPU/model comparisons are
pending. Source/cleanup audits do not establish runtime or BFB results.
