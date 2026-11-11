# Limit turbulent fluxes and compensate the mean field

This draft implements all seven Native `mono_flux_limiter.F90` routines in
source order. The two public interfaces keep all26 and13 required arguments;
there are no Native optional arguments. All2157Native lines and actual JAX
scope were read before `SOURCE_OUTLINE.json` and implementation. Production
source, CMake and shared generators were not edited. Builds and model
comparisons are pending.

Promote the Native-named header/source together when root reaches this live
boundary. Public physical fields and influencing-level index arrays are
column-first LayoutLeft execution-space Views. Integer grid-level values stay
Native one based; only View addresses subtract one. Per-call arrays use the
existing concrete array owners. No mutable scientific state survives a call.

The main limiter retains the Native parallel screening pass over original
neighbor fluxes. Its integer reduction implements the logical ANY gate. Only
when needed does a column kernel traverse flux levels in physical grid order,
using each preceding updated flux. The source eps test controls subsequent
mean replacement, while the screen itself uses the exact strict limit tests.
The sampled flux-limit arrays retain Native assignment paths and unused-value
sentinel. No extra initialized value masks a missing Native assignment.

Influence-window min/max searches, travel-time searches, the spike integral
and same-column scaling remain sequential in Native direction and order.
Independent columns/cells use Kokkos kernels with adjacent columns contiguous.
Scalar locals are private to the kernel that writes them. Grid component
Views/scalars are captured individually; no owning UDT or stats context is
copied into a device lambda. Loop extents are materialized on the host.

The three private matrix routines preserve Native bands-first
`(3,column,level)` scratch. `mfl_xm_lhs` reuses the stable public
`mean_adv::term_ma_zt_lhs` implementation; that owner must be promoted or
provided alongside this draft. Interpolation and the selected tridiagonal
solver already have live public C++ APIs. No new lower Fortran bridge is
needed. One selected constant is pending in the existing generated model-flags
owner: `l_force_descending_solves=false` (`model_flags.F90:194`). All other
constants, Native debug API and sampled stats APIs already exist.

Packing/transposition occurs only at the actual tridiagonal solver API.
Requested descending solves reverse both bands and rows while packing.
After solving, mutated lhs/rhs and the still-unflipped solution are restored
to their Native post-solver layout before checking all-column fatal status.
The final solution reversal occurs only on success. The parent deliberately
copies changed columns' mean solutions before its own fatal check, matching
Native. It does not copy JAX's reordered fatal/flip behavior.

The source-disabled explicit mean-adjustment and constant-thickness branches
are both retained. JAX omits the constant-thickness body and rewrites the
screening/scan; it also guards square roots and divisors. This draft uses
Native formulas and branches without those extra guards. Native references
to historical MKL support remain source comments; the current Native file
contains no separate active MKL calculation branch to implement.

Stats names, blank names, sampling gates, budgets and update order are
unchanged. Four `field/dt` arguments become explicit call-local Views directly
before the existing host stats boundary. Spike warnings stage integer flags
and print in column order after the kernel, keeping both Native CPU messages.
Fatal messages retain their original placement and include preexisting fatal
columns. No device I/O or scalar f2k/device crossing is introduced.

The separate source audit compares seven definition contracts, both public
headers,110 indexed/scalar assignment trees and24Native calls/128actuals.
The actual solver mapping, logical reduction, whole-column scaling/SUM and
four stats payload adapters are reviewed explicitly. It also checks31lambda
bodies for owning UDT captures, scalar first-use and exact internal call case.
This evidence does not establish compilation or numerical correctness.

No deliberate FMA, forced rounding or arithmetic workaround was added for
BFB. The original solver direction choice is a Native option. Integer squares
expand directly with source grouping and source constants retain their kind.
The copied `sigma_w_i_zm` comment is corrected to standard deviation [m/s],
consistent with its square-root assignment; Native science was not changed.

Root owns hookup, obsolete whole-call bridge removal and the maintained
CPU Fortran versus Serial Kokkos zero-tolerance model gate. Nondefault solve
direction, grid direction, mean-advection/spike options, disabled source
branches and fatal paths currently have source review only. GPU execution
and performance are later checkpoints; no runtime or speed claim is made.
