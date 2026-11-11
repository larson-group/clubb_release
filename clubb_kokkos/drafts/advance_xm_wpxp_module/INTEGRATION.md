# Mean/flux advance integration

The public `clubb_kokkos::advance_xm_wpxp` retains all 107 native argument names,
order and required status. Inputs/outputs use the established execution-space
`Real1D`/`Real2D`/`Real3D` types, native grid/implicit-coefficient/resolution types,
caller-owned `err_info_type`, four clipping counters by reference and the live
Fortran `f2k::stats_type`. The existing `advance_clubb_core` caller already uses
these array/type owners; replace its current f2k call only after dependencies
and source review are ready.

Move reviewed source/header to the normal core directory, add the source to the
existing core library, include its header at the live caller and replace the
explicit f2k parent call. This draft makes none of those shared edits. Do not add
the draft folder to production include paths.

## Public lower calls

`lower_contracts.json` contains exact native module/source, argument order,
intent, dimensions and optional status for each dependency. These are real
native public routines; no new physics API or registry is proposed.

| Native module | Required public routine | Arguments |
| --- | --- | --- |
| advance_helper_module | calc_wp3_on_wp2 | 8 |
| mono_flux_limiter | calc_turb_adv_range | 13 |
| diffusion | diffusion_zm_lhs | 10 |
| mean_adv | term_ma_zm_lhs | 7 |
| diffusion | diffusion_zt_lhs | 10 |
| mean_adv | term_ma_zt_lhs | 10 |
| turbulent_adv_pdf | xpyp_term_ta_pdf_rhs | 12 |
| turbulent_adv_pdf | xpyp_term_ta_pdf_lhs | 12 |
| turbulent_adv_pdf | xpyp_term_ta_pdf_lhs_godunov | 8 |
| mono_flux_limiter | monotonic_turbulent_flux_limit | 26 |
| pos_definite_module | pos_definite_adj | 10 |
| fill_holes | fill_holes_vertical_api | 10 |
| clip_explicit | clip_covar | 7 |
| stats_netcdf | stats_update_budget_2d | 3 |

All those lower arguments are required; none is optional. The stats wrapper
should expose the existing generic name `f2k::stats_update_budget`, accepting
`name`, rank-2 `values` and live `stats`, as its existing begin/finalize siblings
do. Add missing bridges through the current generator/manifest owner, then
remove this draft's declaration-only `f2k_dependencies.hpp`. Its declarations
are integration contracts, not implemented bridges or fallback code.

The mean-advection weight arrays retain their native signed third-axis bounds:
`m_above:m_below` or `t_above:t_below`; their C++ owner packs the same storage.
Limiter range outputs and hole-fill boundary indices retain native 1-based
values at the Fortran boundary. Generic View fields/subviews need the existing
pack/unpack owner, not whole-type or global storage.

Existing C++ grid operators, `Skx_func`, error status checks and band solver are
called directly. Stats update/begin/finalize remain on the current Fortran
backend. Of 79 native stats calls, 71 are retained in source order; eight belong
to the excluded sponge blocks. Temporary array-valued budget expressions are
computed explicitly on device, then passed to the existing stats bridge. The
rank-1 reciprocal condition-number diagnostic uses `f2k::Real1D(values)` to
select the existing const execution-View overload; the stats owner stages it.

## Constants and supported case guard

The existing generated `model_flags.hpp` owner must expose these native
compile-time flags before compilation:

- `l_pos_def=false`, `l_clip_turb_adv=false`: model_flags.F90:107–108.
- `l_explicit_turbulent_adv_wpxp=false`: model_flags.F90:122.
- `l_force_descending_solves=false`: model_flags.F90:194.

The draft includes the corresponding `clip_explicit.hpp` constant owner, which
is not currently present in production. Establish its native named constants
through the existing constant-owner approach when that module is integrated:

| Native constant | Value | Source line |
| --- | --- | --- |
| clip_wprtp | 8 | clip_explicit.F90:25 |
| clip_wpthlp | 9 | clip_explicit.F90:26 |
| clip_upwp | 10 | clip_explicit.F90:27 |
| clip_vpwp | 11 | clip_explicit.F90:28 |
| clip_wpsclrp | 13 | clip_explicit.F90:30 |
| wprtp_cl_max, wpthlp_cl_max, upwp_cl_max, vpwp_cl_max | 3 each | clip_explicit.F90:37–40 |

No shadow constants were added to this physics draft.

Actual JAX module lines 23–24 explicitly exclude sponge damping because its
case initialization rejects sponge-enabled cases. Native global settings for
rtm, thlm and uv sponge damping are not C++ arguments to this routine. Root must
retain that supported-case rejection at the established initialization/bridge
boundary before hooking up this parent; comments alone do not enforce it. The
source marks all three exact native boundaries, immediately before prognostic
wind nudging. No sponge state or replacement helper is invented here.

## Necessary language and execution adaptations

- Private templates accept actual field/scalar subviews without reallocating
  their data or widening the public API. Definitions remain in native source
  order; forward declarations satisfy C++ lookup.
- Native loop variables remain 1-based within routine calculations; View
  indexing subtracts one. The existing zero-based tunable-parameter owner is
  used directly for its parameter axis. Array assignments expand into explicit
  cell loops; scalar/field loop order around lower calls remains native.
- Independent loops execute on Kokkos. Scratch arrays, including rank-4 scalar
  LHS coefficients, use LayoutLeft with no mutable global state. Native physical
  fields remain `(column,level[,scalar])`; private matrix arrays retain native
  `(band,column,row)` indexing for source comparison.
- The established matrix wrapper expects `(column,band,row)`. The private
  `xm_wpxp_solve` transposes those two axes at that one boundary, also applying
  native descending row/band reversals when requested. RHS reversal snapshots
  are allocated only when reversal is requested. Solution reversal uses a
  separate result to preserve native whole-array assignment alias semantics.
  Mutated LHS/RHS and optional condition estimates return to caller storage.
- Genuine optionals use pointers: `Real1D* rcond=nullptr` and
  `int* wpxp_cl_num=nullptr`. No new optional fields were introduced.
- Array ANY conditions use the existing status owner or a small integer Kokkos
  reduction. Clipping correction uses an ordered per-column mask scan with an
  integer column count; its undefined source entries are documented separately.
- Native host I/O stages fields through `array_types.hpp`. Diagnostic arrays
  print in native first-index-fastest order with max_digits10 precision. Fixed-
  width error headers retain the established `(character,column)` shape. C7
  validation and hole warnings snapshot their input once for the host loop.
  Printing, bridge calls and host reads never occur inside kernels.
- The explicit scalar TA rank mismatch uses a documented unmanaged rank-2
  alias of the first native contiguous rank-3 LHS elements, preserving native
  sequence association rather than silently assigning the JAX RHS repair.

Small adapters for array ANY and host diagnostic formatting contain no
field-specific physics or hidden call lists. Matrix and native sequence
association changes are data-layout boundaries. No speedup is claimed before
measurement. Native source-definedness choices and BFB policy are in
`BFB_NOTES.md`.

The `wpxp_term_tp_lhs::lhs_tp` comment corrects a copied Native unit label:
`wp2 * invrs_dzm` has units `[m/s^2]`, rather than `[1/s]`. The source
assignments and executable tokens are unchanged by the formatting cleanup.

## Validation handoff

No build, compilation, install or numerical test was run for this draft.
Integrate a reviewed portion with its declared lower dependencies, then use the
maintained short `multi_col_stats` runner with CPU Fortran/CPU Kokkos zero
tolerances. Broader duration/stats checks belong to larger integration gates;
GPU comparisons are occasional milestones. Record any observed answer change
and any deliberate arithmetic adjustment in the existing BFB adaptation owner.
