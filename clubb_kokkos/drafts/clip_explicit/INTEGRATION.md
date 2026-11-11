# Integration contracts

Read [DRAFT_STATUS.md](DRAFT_STATUS.md) first. Promote reviewed portions into
the Native-named production owner only after the active CPU comparison passes.
Root owns shared source registration and obsolete whole-routine bridge removal.

## Public owner and interfaces

[`clip_explicit.hpp`](clip_explicit.hpp) defines all 19 Native public integer
constants, retaining the gaps in the selector numbering. Native `clip_scalar=7`
remains commented out. Promote these definitions as the shared clipping owner;
do not duplicate selected constants in each moment parent.

All six public interfaces retain Native argument names, order, intents and
optional status:

| Routine | Arguments | Mutation |
| --- | ---: | --- |
| `clip_covars_denom` | 24 required | Stats, four counters and covariance fields |
| `clip_covar` | 7 required | Covariance; Out change field |
| `clip_variance` | 8 required, 1 optional | Stats and variance |
| `clip_skewness` | 10 required | Stats and third moment |
| `clip_skewness_core` | 8 required | Third moment |
| `clip_rcm` | 5 required | Cloud-water field |

The optional `threshold_hi` is a scalar `const real_t*` with null meaning absent.
No optional arrays or extra error-state arguments are introduced. Fields are
column-first LayoutLeft execution-space Views and retain caller allocations.
Grid metadata remains Native one based. Scratch is local to each call.

The public generic `clip_covar` and `clip_variance` overloads are justified by
actual scalar-section callers in the stable moment drafts. They reuse the
existing array owner to adapt memory/layout once at entry and copy back only
Out/InOut fields. Compatible execution-space Views alias existing allocations.
The covariance change field is a fully written Native Out field, so its entry
adaptation does not copy old values. The other four public interfaces use the
existing concrete execution-space View types.

## Dependencies and statistics

The production constants owner already supplies `max_mag_correlation`,
`max_mag_correlation_flux` and `zero_threshold`. No additional model constants
are required by this module.

Reuse the stable
[`advance_helper_module.hpp`](../advance_helper_module/advance_helper_module.hpp)
owner for `smooth_heaviside_peskin(nz,ngrdcol,input,smth_range)`. The function
returns an owned rank-two execution-space View. Its result is copied into the
Native-named local `H_zagl`; do not introduce a private or duplicate lower
implementation here.

Statistics remains Native. All 20 budget call sites preserve their original
order, names and gates. A Native blank budget name is still passed when the
source passes it; JAX's extra nonblank-name guard is not copied. Array actual
arguments such as `xp2/dt` become call-local Views immediately before the stats
API boundary. No stats object enters a device kernel.

The one missing public lower definition is
`stats_update_budget(name,values,stats)` for Native
`stats_netcdf.F90:stats_update_budget_2d`, with all three arguments required and
a rank-two In field. [f2k_dependencies.hpp](f2k_dependencies.hpp) contains only
a pending declaration; [lower_contracts.json](lower_contracts.json) records the
actual Native contract. Extend the existing transient stats/bridge owner, then
remove this staging declaration. There is no raw/private binding or fallback.

## Diagnostic and termination boundaries

`clip_variance` retains the optional invalid-threshold guard absent in JAX.
An integer device reduction identifies `threshold_lo > threshold_hi` over
numeric levels `1..nzm-1`. Only failure stages the threshold field to the host;
the host prints the Native values at the first failing source-order column and
level and executes `std::exit(EXIT_SUCCESS)`. Native unqualified `STOP` is normal
process termination, not an `err_info` update or returning exception. No budget
or clipping pass precedes this guard. Parallel Native GPU diagnostic ordering
is not specified; the host failure report follows the Native serial loop order.

At debug level 3, variance warning passes use packed host mirrors and the
original compare/write/assignment order. Lower warnings print before assigning
the minimum; upper warnings print the already-assigned maximum, exactly as the
source does. Ordinary variance clipping runs device kernels without host
staging. `clip_rcm` stages its warning inputs only under its Native debug gate;
all warnings precede its actual device update. Text payloads, values and indices
are preserved; C++ numeric formatting uses full precision rather than Fortran
list-directed formatting.

## Source-definedness and branch details

The full JAX file was read for scope; Native remains authoritative:

- Covariance clipping uses strict upper/lower comparisons. Neither boundary
  covariance is altered, and both boundary changes are zero. If a negative
  variance product makes the bound NaN, neither comparison is true: covariance
  remains unchanged and change is zero. Do not replace this with a NaN-producing
  unconditional clip or add an invented variance floor.
- All four covariance counters increment regardless of prediction flags.
  Perturbed wind calls reuse the same change scratch after the base calls,
  preserving the source overwrite order. Scalar clipping has no new budget.
- Minimum variance clipping follows physical grid direction and excludes the
  physical upper boundary. Optional maximum validation and clipping retain
  the different numeric `1..nzm-1` bounds in either grid orientation. No attempt
  is made to repair this Native asymmetry silently.
- Both default and smoothed height-dependent skewness limits are included,
  with the `0.0021` factor, 100 m altitude boundary and absolute `wp3=100` cap.
  Negative variance/limit square roots are not guarded beyond Native behavior.
- `clip_rcm` subtracts the floating kind's `EPSILON`, not a spacing depending
  on `rtm`. Its source limitation for nonpositive total water remains explicit;
  no additional water-conservation repair is introduced.
