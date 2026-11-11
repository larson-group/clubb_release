# Integration contracts

Read [DRAFT_STATUS.md](DRAFT_STATUS.md) before promoting this isolated draft.
Root owns shared source registration, lower bridges and constants, removal of
the obsolete whole-routine bridge, and CPU bit-for-bit verification.

## Public interface and ownership

[`advance_windm_edsclrm_module.hpp`](advance_windm_edsclrm_module.hpp) retains
all 52 Native arguments, in source order, with no new optionals. Scalar clipping
counters are references; field outputs keep the caller's allocation. Physical
fields use column-first LayoutLeft execution-space Views. Grid/scalar metadata
retains Native one-based values; View addresses subtract one. The existing
generated tunable parameter indices are already zero based.

All scratch is per call and uses `WithoutInitializing`. The two explicit-copy
validators execute before device scratch allocations. Each forcing scan retains
one serial vertical scan per independent column, so status updates do not race.
Both validators mark only affected columns, check the complete existing status
array and immediately return on any fatal error, including a preexisting fatal
column. Status is not reset at entry.

Private helper View templates accept the source rank-preserving scalar sections.
The five private routines and the main control flow are implemented here. The
only target-only helper formats complete arrays, in first-index-fastest order,
at the existing host diagnostic boundary. It reuses `to_host_space` to pack
strided Views. No host conversion occurs inside a physics kernel.

## Dependencies

Reuse existing production grid, matrix solver, UDT/status, tunable and stats
owners. Promote the stable
[`advance_helper_module.hpp`](../advance_helper_module/advance_helper_module.hpp)
owner for its Native `calc_xpwp` generic and `pvertinterp`; do not add duplicate
flux/interpolation helpers. Both retain caller-owned Out Views. Scalar-field
sections use the existing generic View adaptation of `calc_xpwp`.

[`f2k_dependencies.hpp`](f2k_dependencies.hpp) contains declaration-only pending
public lower calls, with complete Native contracts in
[lower_contracts.json](lower_contracts.json). These are not implementations or
fallbacks. Add definitions through the existing transient `src/f2k` owner, or
promote the matching C++ owner, then remove the staging declarations.

| Public lower call | Native owner | Required arguments | Relevant result |
| --- | --- | ---: | --- |
| `diffusion_zt_lhs` | `diffusion.F90` | 10 | Out `lhs(3,ngrdcol,nzt)` |
| `term_ma_zt_lhs` | `mean_adv.F90` | 10 | Out `lhs_ma(3,ngrdcol,nzt)` |
| `clip_covar` | `clip_explicit.F90` | 7 | InOut covariance; Out change |
| `fill_holes_vertical_api` | `fill_holes.F90` | 10 | InOut rank-two field |
| `stats_update_budget` | `stats_netcdf.F90:stats_update_budget_2d` | 3 | Native stats context |

All listed arguments are required. In particular, the lower LHS arrays retain
Native bands-first order. Reuse these same contracts already pending for the
other timestep drafts rather than generating separate bridge families.

Two existing constant owners need their selected Native definitions:

- `model_flags.hpp`: `l_force_descending_solves = false`, from Native
  `model_flags.F90:194`.
- `clip_explicit.hpp`: `clip_upwp = 10`, `clip_vpwp = 11`,
  `upwp_cl_max = 3` and `vpwp_cl_max = 3`, from Native
  `clip_explicit.F90:27–40`. This header is not currently present in production.

Keep these in their Native-named/generated owner, with no shadow physics-file
constants or raw/private bindings. The three `windm_edsclrm_*` selectors are
Native private module constants and remain local here.

## Matrix and statistics boundaries

Native LHS scratch is `(band,column,row)`; the established C++ tridiagonal API
uses `(column,band,row)`. `windm_edsclrm_solve` transposes only at that actual
solver call boundary. Its packing also applies the original optional band/row
reversal. Only `nrhs` RHS/solution vectors are touched, preserving unused storage
when the main caller has `max(2,edsclr_dim)` vectors.

After the solver call, restore its mutated LHS/RHS and unflipped solution into
the Native caller allocations before the Native debug/fatal check. A fatal
return therefore retains the source side effects. Only a non-returning path
reaches the original final solution-direction reversal. Optional `rcond`
dispatch and reciprocal-condition statistics remain in Native order. The
existing solver's supported method scope still applies.

Statistics remains Native. All 24 supported field-specific update/budget call
sites retain their original order and gates, including the source clipping
budget guards inside the `!l_predict_upwp_vpwp` branch. Arithmetic actual
arguments such as `field/dt` and `1/rcond` are evaluated into call-local Views
immediately before the public stats call. Rank-one reciprocal-condition output
explicitly selects the existing `f2k::Real1D` overload. OpenACC host-only update
blocks become comments: the stats API stages its own inputs.

## Supported branches and source-definedness

The draft includes both grid orientations, standalone/host forcing behavior,
LMM stepping, nudging, predicted versus diagnosed momentum fluxes, perturbed
winds, arbitrary eddy scalars, explicit thlm/rtm copy diffusion, clipping,
hole filling, sampled statistics and Native diagnostics/fatal returns.

The UV sponge block is the one scope exclusion. JAX rejects it in
`clubb_case_initalization._check_unsupported_features`; the current C++
[`clubb_driver.cpp`](../../src/clubb_driver.cpp) rejects all sponge settings
before allocation. Retain that gate for every new caller; the core draft has
no sponge state argument. Its exact omitted block is marked at the source
position. Supporting sponge cases requires porting that original block and
its owner rather than silently allowing the setting.

C++ follows Native rather than JAX for these distinctions:

- Validators and wind-solve fatal errors return immediately. The eddy-scalar
  solve prints its fatal message, continues to the full source field dump,
  then returns at the original final debug check.
- `wind_speed` and `u_star_sqd` are unwritten when momentum fluxes are predicted.
  The scalar LHS receives their handles with implicit surface flux disabled and
  never reads their elements. JAX initializes these otherwise-unused fields.
- `um_old`/`vm_old` are written only in the diagnosed-wind LMM branch. A later
  fatal scalar dump with predicted winds and LMM can print their unwritten
  values, as in Native. This is not a defined-output BFB claim or an attempted
  repair. The draft does not add JAX's zero initialization.
- Scalar flux first uses `Km_zm_p_nu10`, while its later assignment uses
  `Kmh_zm` and replaces the interior value. The draft retains both the distinct
  coefficients and replacement, instead of combining two halves speculatively.
