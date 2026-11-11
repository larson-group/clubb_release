# Integration contracts

Read [DRAFT_STATUS.md](DRAFT_STATUS.md) first. Promote reviewed portions into the
Native-named production owners only after the active CPU comparison passes.
Root owns CMake registration, public lower implementations, missing constants,
upstream scope guards and obsolete whole-routine bridge removal.

## Public and private interfaces

`advance_xp2_xpyp_module.hpp/.cpp` retains all 66 required timestep arguments and
all 18 required `update_xp2_mc` arguments. Physical fields use column-first
LayoutLeft Views, and every Out/InOut field keeps the caller's allocation.
Private procedures use rank-preserving View templates so Native scalar slices
can be passed without target-only physics wrappers. All sixteen Native routines
remain in source order; only the two public routines leave the unnamed
namespace. There are no new optional arguments.

Native scalar/grid category metadata remain one based. View addresses subtract
one; generated `parameter_indices.hpp` indices are already zero based and are
not shifted a second time. Three legacy declarations omit INTENT:
`term_dp1_lhs::lhs` is OUT, `term_dp1_rhs::threshold` is IN and
`term_pr2::rhs_pr2` is OUT. Their complete Native bodies establish these contracts;
only the draft outline/parser records overrides.

## Explicit lower boundaries

[f2k_dependencies.hpp](f2k_dependencies.hpp) contains declaration-only pending
public lower calls, with exact Native contracts in
[lower_contracts.json](lower_contracts.json). It implements no bridge or fallback.
Promote an available C++ lower routine or extend the existing temporary bridge
owner, then remove these staging declarations. No new private binding API or
global storage is needed.

Required lower names are `calc_wp3_on_wp2`, `diffusion_zm_lhs`, `term_ma_zm_lhs`,
`xpyp_term_ta_pdf_lhs`, `xpyp_term_ta_pdf_rhs`, both corresponding Godunov routines,
`fill_holes_vertical_api`, `clip_covar` and `clip_variance`. The complete
advance-helper draft already contains `calc_wp3_on_wp2`; promote/use that owner
where appropriate rather than adding a permanent duplicate.

`clip_variance` has the genuine Native optional scalar `threshold_hi`; none of
this module's calls supplies it. The pending declaration retains absent/present
semantics with a nullable scalar pointer. The other listed science arguments
are required. Root must preserve Native bands-first matrix ranks when extending
these lower APIs or replacing them with C++ ports.

Existing grid APIs, `matrix_solver_wrapper`, status/UDT owners and stats updates
are reused directly. No whole `advance_xp2_xpyp` or private routine is delegated
to Fortran.

## Native statistics backend

All field-specific begin/update/finalize calls and conditional registrations
remain in Native order. Array expressions such as `field/dt`, reciprocal
condition numbers and scaled LHS actual arguments are materialized into local
Views immediately before their calls. These are language array-expression
buffers, not shared stats state. Rank-one condition output explicitly selects
the existing `f2k::Real1D` overload; rank-two updates use the public generic
stats boundary.

Two missing stats capabilities need the existing owner:

- Public `stats_update_budget(name,values,stats)` for the Native 2D routine.
- `stats_finalize_budget(name,values,stats,l_count_sample)` with the true Native
  optional logical. Three clipping calls explicitly supply false so an
  additional clipping correction does not increment the budget sample count.
  Existing three-argument calls must keep their true default. The draft declares
  a four-argument explicit-value overload only; root may implement one unified
  optional-capable public interface when promoting it.

## Matrix storage and ordered work

Native variance matrices use `(band,column,row)`; the established C++ solver
uses `(column,band,row)`. Private scratch deliberately retains Native order to
match every lower diffusion/advection/stats contract. `xp2_xpyp_solve` visibly
transposes only at the existing solver boundary. It combines the Native optional
row/band reversal with that packing. RHS and solution remain
`(column,row,RHS)`; rank-two single-RHS callers use the same Native sequence
association through explicit rank dispatch. Only `nrhs` requested vectors are
read/written, preserving untouched extra storage in the two-RHS wind buffers.

The chosen solver mutates the Native InOut lhs/rhs allocations before the Native
fatal diagnostic. Its unflipped result is also restored before that diagnostic;
only a successful return reaches the original final answer reversal. Optional
`rcond` dispatch, shared-LHS condition names and all-column error tests remain
unchanged. The existing solver's supported method scope still applies.

The single-LHS eligibility test reduces integer disagreement hits and retains
Native `abs(C2a-C2b) > abs(C2a+C2b)/2*eps` expressions. No floating reduction is
introduced. `update_xp2_mc` diagnoses precipitation fraction independently by
column while retaining its serial physical-downward level recurrence in both
grid orientations. Other independent cells map to Kokkos kernels with adjacent
columns fastest where the Native matrix layout exception does not apply.

## Constants and exclusions

Use the existing generated model/parameter constants. `l_force_descending_solves`
is still missing from the selected model-flags owner, as already required by
other timestep drafts. The existing Native `clip_explicit` constant owner also
needs its C++ header, including `clip_sclrp2=14`, `clip_sclrprtp=15` and
`clip_sclrpthlp=16`. Do not add shadow constants to the physics file. No new
physical coefficient is needed.

Sponge-enabled up2/vp2 cases are excluded with JAX scope. The single Native
sponge block and its four sampled budget calls are marked at their exact source
location; initialization must reject those settings before this draft runs.
This is an explicit unsupported option, not a successful damping port.

No newly defined scientific fallback repairs Native-unwritten states. The
initialized-false `l_clip_large_neg_mc` branch remains in source order, including
its existing diagnostic temporary lifetime. Unknown PDF selections retain the
Native branch surface; this draft does not invent implicit turbulent-advection
coefficients for source-unhandled selections. These facts are not numerical
validation for dormant or source-undefined inputs.
