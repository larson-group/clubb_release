# PDF utility draft

[pdf_utilities.hpp](pdf_utilities.hpp) implements all 19 public routines from
[Native pdf_utilities.F90](../../../src/CLUBB_core/pdf_utilities.F90), in Native
definition order. It is source-audited and formatted, but has **not been
compiled, integrated, executed or checked for BFB**. No production files,
CMake inputs or generated bridge manifests were changed.

The header contains 15 device-callable scalar elemental definitions and four
array procedures. Header-only definitions keep scalar formulas available to
GPU callers and preserve the complete Native routine order. Four distinct
`_dp` routines are included even though JAX omits those duplicates: the task
covers the complete Native module, and `mean_L2N_dp` deliberately lacks the
core version's `tiny` guard. The double correlation routines also retain
Native's conversion of the core-kind correlation limit to double.

## Integration

1. Promote the PDF header into `clubb_kokkos/src/CLUBB_core/`. Include it in the
   existing PDF caller and replace its reached `f2k::` calls with the same
   Native-named C++ interfaces. No new C++ compilation unit is needed.
2. Merge the exact scalar `smooth_max_scalars` definition and scalar
   `smooth_max` generic from [advance_helper_module.hpp](advance_helper_module.hpp)
   into that module's eventual permanent owner, in Native procedure order.
   This partial helper header covers only Native lines 1723–1761; other
   `advance_helper_module` procedures remain outside this draft. Root owns
   promotion and dependency deduplication.
3. Remove the obsolete whole-routine bridge entries for
   `calc_comp_corrs_binormal` and `smooth_corr_quotient` when their final C++
   callers are replaced. Remove the transient PDF array adapter and its CMake
   source entry when no `f2k::compute_variance_binormal`, `calc_corr_chi_x` or
   `calc_corr_eta_x` callers remain. Keep unrelated bridge entries intact.
4. Compile and use the maintained short multicolumn model comparison for
   acceptance and diagnosis. Native backend calls made by still-unported
   routines remain Native; this draft changes none of them.

The existing rank-two elemental calls have the same argument order and return
new execution-space `Real2D` values. Their three overloads cover only
`compute_variance_binormal`, `calc_corr_chi_x` and `calc_corr_eta_x`; scalar
forms remain directly callable inside GPU kernels. Inputs may be strided but
must be conformable and accessible from the default execution space. The four
Native array procedures retain explicit dimensions and caller-owned outputs.
No mutable global storage, derived-type registry, stats interface or scalar
Fortran crossing is introduced.

The only lower dependency is `smooth_max_scalars`; it is included exactly in
the partial helper header. Required constants already exist in
`constants_clubb.hpp`: `one`, `one_dp`, `zero`, `zero_dp`, `two`, `one_half`,
`eps`, `max_mag_correlation`, `min_max_smth_mag`, `chi_tol`, `eta_tol`,
`rt_tol` and `thl_tol`. Native `tiny` maps to the minimum normalized value of
`real_t` from `std::numeric_limits`.

## Review and numerical status

[source_outline.json](source_outline.json) records all contracts and branches.
[source_audit.json](source_audit.json) records the separate source review:
all 19 argument lists and definition order, all assignment expression trees
and branch predicates, plus the scalar helper's expression tree. Every
integer square is checked with its surrounding arithmetic grouping intact.
Manual review covers scientific comments/citations, precision differences,
strict tolerance comparisons, ordered clipping, each `calc_xp2` branch,
sequential smoothing calls, array ownership and the three View adaptations.
The module has no Native diagnostics or stats updates to reproduce.

No explicit FMA, forced floating-point association, added floor, `log1p`,
`expm1` or compiler-specific math workaround was added. Native integer
squares are parenthesized products. BFB remains pending; any demonstrated
adjustment after integration belongs in the permanent BFB adaptation log
with its Native compiler evidence and maintained model comparison result.

The preparation and source-audit scripts are draft artifacts. They prepare
or compare source only; they do not execute physics or establish numerical
equivalence and should not be copied into the production module.
