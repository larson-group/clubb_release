# Saturation integration

Move the reviewed `saturation.hpp` into the existing core include directory when
ready to replace its live Fortran boundary. Scalar definitions remain in the
header because device parcel callers need visible definitions. There is no
separate compiled translation unit, new mutable state or Fortran dependency.
Do not add this draft directory to production include paths.

## Public contracts

All APIs are in `clubb_kokkos`. Scalar inputs/results use the existing `real_t`.
Pressure is in Pa, temperature in K, vapor pressure in Pa and mixing ratios in
kg/kg. Formula selectors are native integer constants.

```cpp
// All four scalar APIs are KOKKOS_INLINE_FUNCTION.
real_t sat_mixrat_liq_api(real_t p_in_Pa, real_t T_in_K, int saturation_formula);
real_t sat_mixrat_ice(real_t p_in_Pa, real_t T_in_K, int saturation_formula);
void sat_vapor_press_liq(real_t T_in_K, int saturation_formula, real_t& esat);
real_t rcm_sat_adj(real_t thlm, real_t rtm, real_t p_in_Pa, real_t exner,
                  int saturation_formula);

// Pressure and Temperature are rank-2 Kokkos Views indexed (column, level).
auto sat_mixrat_liq_api(int nz, int ngrdcol, const grid& gr,
                       Pressure p_in_Pa, Temperature T_in_K,
                       int saturation_formula,
                       std::optional<int> start_index_in = std::nullopt);
auto sat_mixrat_ice(int nz, int ngrdcol, Pressure p_in_Pa, Temperature T_in_K,
                    int saturation_formula);
```

Array calls return owned, LayoutLeft real Views of shape `(ngrdcol,nz)` in
`Pressure::memory_space`. Existing `array_types.hpp` owns input staging and
output transfer. No arbitrary host/device deep copy is duplicated here. Passing
an integer start index converts to the optional; omission retains the native
full-grid behavior. Ascending grids select `start_index_in:nz`; descending grids
select `1:start_index_in`. Indices at this public boundary remain native 1-based.
Outside the selected liquid range, output is unwritten/undefined in native and
in this draft. An empty native range launches zero cells.

The mixing-length draft's required three-scalar parcel call and seven-argument
array call have these exact positional contracts. Public PDF closure uses the
same liquid/ice array generic names as the existing bridge. Call-site and
bridge-removal changes belong to root's integration, not this isolated draft.

## Existing owners and missing constants

Existing headers provide `real_t`, arrays/staging, `grid`, and constants
`ep`, `T_freeze_K`, `Lv`, `Cp`, `zero_threshold`. The installed Kokkos 5.2.2
mathematical-functions header supplies device `exp`, `log`, `log10`, `pow` and
`round`; `Kokkos::abort` provides the device fatal path. No lower physics or
complex-root bridge is needed.

Root must expose these existing native constants through the established
`model_flags.hpp` generator owner before compilation:

| Native name | Native value | Source line |
| --- | --- | --- |
| `saturation_bolton` | 1 | model_flags.F90:160 |
| `saturation_gfdl` | 2 | model_flags.F90:161 |
| `saturation_flatau` | 3 | model_flags.F90:162 |
| `saturation_lookup` | 4 | model_flags.F90:163 |
| `I_sat_sphum` | false | model_flags.F90:185 |

This draft references that owner directly; it does not create shadow constants
or a silent fallback. Include the real production header during integration.

## Necessary adaptations and source behavior

- The immutable native 188:343 lookup table is a private device getter with 156
  literal switch returns. It preserves each source token and needs no mutable
  global initialization, per-call transfer or 156-entry per-thread array.
  Code size and lookup performance remain unmeasured.
- The identical native liquid scalar/2D vapor-pressure blocks share one device
  cell body. The native specific-humidity conversions remain separate because
  their actual denominators differ. Native routine names/order are retained;
  small overloads implement Fortran's two public generic interfaces.
- Scalar ice lifts the native 2D cell calculation directly, replacing its native
  allocation-free 1x1 local-array forwarding path. This permits a scalar device
  call without launching/staging a Kokkos array operation inside each parcel.
  Array ice calls that same scalar body.
- Each array routine fuses independent vapor-pressure and mixing-ratio loops.
  Cells are flattened with adjacent columns first; no vertical dependence or
  reduction exists in this source. Fusion has no measured speedup yet.
- Native ANINT followed by INT is `Kokkos::round` followed by `static_cast<int>`;
  both round halfway away from zero on ordinary finite representable inputs.
  Native index clipping is retained after the conversion.
- Native `error stop` in the scalar adjustment is device `Kokkos::abort` with
  the original message. The exact native iteration/exit/order and unusual final
  loop-counter equality are retained.

`BFB_NOTES.md` records arithmetic choices and inherited source distinctions.
There are no stats calls in this native module.

## JAX comparison and runtime handoff

Actual JAX source deliberately omits lookup, supports only false `I_sat_sphum`
and Flatau ice, and replaces native adjustment iteration with a 64-step scan.
See saturation.py:10–18, 162–166, 169–242 and 245–299. This draft retains all
native selectors and both humidity formulas, without JAX gradient-specific
changes, denominator guards or new clipping.

No compilation or numerical tests were run for this draft. Integrate the next
required boundary and use the maintained short multi_col_stats whole-model
comparison. Any later forced FMA or math-library adaptation requires observed
before/after evidence and a corresponding entry in `BFB_ADAPTATIONS.md`.
