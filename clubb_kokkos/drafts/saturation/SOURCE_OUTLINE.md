# Saturation source outline

All 1185 lines of `src/CLUBB_core/saturation.F90` were read before implementation.
The native public APIs are liquid mixing ratio, ice mixing ratio, scalar vapor
pressure and scalar saturation adjustment. The scalar forms need device calls
from mixing-length parcel trajectories.

| Native routine, in source order | Lines | Arguments |
| --- | --- | --- |
| `sat_mixrat_liq_k` | 88–280 | p_in_Pa, T_in_K, saturation_formula |
| `sat_mixrat_liq_2D` | 284–549 | nz, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula, optional start_index_in |
| `sat_vapor_press_liq` | 553–621 | T_in_K, saturation_formula, output esat |
| private `sat_vapor_press_liq_flatau` | 624–723 | T_in_K, output esat |
| private `sat_vapor_press_liq_bolton` | 727–760 | T_in_K, output esat |
| private `sat_vapor_press_liq_gfdl` | 765–814 | T_in_K, output esat |
| private `sat_vapor_press_liq_lookup` | 818–857 | T_in_K |
| `sat_mixrat_ice_k` | 861–908 | p_in_Pa, T_in_K, saturation_formula |
| `sat_mixrat_ice_2D` | 912–1086 | nz, ngrdcol, p_in_Pa, T_in_K, saturation_formula |
| `rcm_sat_adj` | 1089–1183 | thlm, rtm, p_in_Pa, exner, saturation_formula |

Preserve native selectors 1 Bolton, 2 GFDL, 3 Flatau and 4 lookup. Invalid
selectors use the -99999.999 vapor-pressure sentinel; they do not throw new
exceptions. Ice has no lookup case and uses the sentinel for selector 4.
`I_sat_sphum` is a compile-time native flag, currently false. Retain both
branches and the actual scalar/2D liquid denominator difference.

Liquid Flatau uses the factored eighth-order polynomial, -85 deg_C floor,
Table 3/4 references and issue 834 explanation. Bolton uses the exponential.
Liquid mixing-ratio GFDL has a 173.15 K floor; public vapor-pressure GFDL instead
has a 203.15 K floor. Lookup has 156 immutable values for 188–343 K, ANINT ties
away from zero, then a clipped integer index. Scalar/2D liquid have identical
vapor-pressure cell calculations; their specific-humidity conversions differ.

Liquid 2D optional `start_index_in` retains native 1-based indexing and the
direction-dependent selected range. Output outside that range is unwritten in
native. Array APIs return owned storage in the pressure field's memory space.
All cells are independent: no vertical recurrences, stats or persistent state.

Ice scalar forwards through native 1x1 local arrays. The device form must lift
the same cell body without an array API call inside each parcel. Ice 2D has
Bolton log/exponential, Flatau Horner with a -90 deg_C floor, GFDL with a
173.15 K floor, then pressure guard and the two humidity formulas.

`rcm_sat_adj` has tolerance 0.001 K, itermax 1000000, ordered answer/convergence/
bracket updates, first-iteration upper-bound override, and final equality test
`iteration==itermax`. Preserve native counted-loop semantics, including natural
exhaustion. Device error stop maps to `Kokkos::abort` with the original message.

Actual JAX source omits lookup, uses false `I_sat_sphum`, supports Flatau ice
only and replaces adjustment iteration with a 64-step scan. This draft covers
all native selectors/branches without JAX gradients, extra guards or new floors.
Math dependencies are Kokkos exp/log/log10/pow/round and abort; installed Kokkos
5.2.2 headers confirm round exists. No lower physics bridge is needed. The
existing generated model_flags owner must expose the native selectors and
humidity flag during root integration; this task does not edit shared owners.

Adaptations to audit: immutable device table getter; shared liquid vapor-
pressure cell body; scalar ice extraction; fused array vapor-pressure/conversion
loops; memory staging; device abort and MIN/MAX translation. No speculative
FMA or BFB workaround is planned. All runtime evidence remains pending.
