This isolated draft ports all 18 procedures in `src/CLUBB_core/setup_clubb_pdf_params.F90`.
It is ready for independent source review and incremental integration. It has not
been compiled or compared with a model run. No CPU, GPU or BFB claim is made.

`SOURCE_OUTLINE.json` records the complete Native and actual JAX reads before
implementation. `SOURCE_AUDIT.json` records the separate owner source review;
`CLEANUP_AUDIT.json` records the formatting and comment pass. No independent
full-file review has been completed for this draft.

The six public interfaces and twelve private helpers preserve Native procedure,
argument, field, branch and call order. The parent owns scratch for precipitation
fractions, component moments, transformations and matrix factors. Stats uses the
existing Native collector. Host matrix diagnostics retain their existing APIs.

To integrate, copy the source/header into the corresponding core directory and
satisfy these existing-owner dependencies before removing the whole-procedure
bridge:

- Extend `array_types.hpp` with `Real4D` and rank-4 allocation branches in its
  existing host/execution conversion functions. Correlations and Cholesky factors
  use `(column,level,pdf_index,pdf_index)` with `LayoutLeft`. This draft deliberately
  does not create a competing array owner or duplicate rank-4 conversion helpers.
- Promote the stable `matrix_operations`, `diagnose_correlations_module`,
  `corr_varnce_module` and `index_mapping` drafts. The metadata procedure fragments
  must merge into the existing generated `hm_metadata_type` owner. The matrix
  owner also requires its recorded LAPACK shim extensions.
- Promote `clip_explicit` for `clip_covar` and `clip_wphydrometp`. Existing live
  grid, PDF utilities and advance-helper APIs supply the remaining scientific
  calls. All referenced physical and tunable constants already exist in their
  live headers; this file adds only its Native private boolean parameter.
- Implement the two pending transient array contracts declared in
  `f2k_dependencies.hpp` in the existing f2k owner, or promote their actual Native
  owners and replace those two call sites. `precip_fraction` has 22 arguments,
  including its Native error and stats objects. `Nc_in_cloud_to_Ncnm` is an
  elemental Native function with ten real inputs; its temporary shape adapter
  adds `nzt,ngrdcol`, evaluates a conformable whole-array assignment in Fortran,
  and returns an owned HostSpace rank-2 View. There is no per-cell bridge call.
  No bridge implementation is claimed here.

The public array adapters and ten private section adapters reuse
`to_execution_space`/`from_execution_space`. They preserve caller allocations and
restore every Native Out/InOut array. In particular, the per-cell rank-4 matrix
sections are strided and must pass through the matrix owner's existing adapters.
The optional hydrometeor output remains in its original argument position as a
nullable pointer to the existing `hydromet_pdf_parameter_collection`. The caller
must allocate that collection and the precipitation fields through their existing
owners; this port does not invent an allocation registry. Only active hydrometeor
slots are assigned, as in Native.

Device kernels capture individual component Views or scalar metadata. The two
Native elemental helpers are device-callable. Independent array assignments and
explicit cells become column-first kernels; source species and triangular
assignment stages retain their order. The Native per-column elemental `rcm_pdf`
assignment becomes one independent column/level kernel. Metadata lookups and
character handling stay on the host with call-local mirrors. The small TRIM
adapter removes trailing blanks only. There is no mutable global storage, device
host callback, scientific wrapper around the parent or speculative FMA change.

The following Native/JAX boundaries are deliberate and need to remain reviewable:

- Native lines 573–574 pass `sigma_w_1` for both formal component standard
  deviations in `compute_mean_stdev`; this draft preserves that actual call.
  JAX passes `sigma_w_2` as the second actual. Correcting Native is separate work.
- Native normal transforms use negative `HUGE(real_t)` for absent means. Normal
  mean stats clamp those values to negative default-real `HUGE(0.0)`, widened to
  `real_t`, before output. The draft preserves both kinds and all stats guards;
  JAX's differentiability floor is not imported.
- The nonfixed eta–w branch in `comp_corr_norm` rewrites the w–chi slots in
  Native. The draft preserves this, as JAX does.
- `denorm_transform_corr` writes the diagonal and lower triangle only. Its upper
  entries are undefined in Native. JAX symmetrizes them; this draft does not.
- Native computes true-space correlations only when the optional hydrometeor
  output is present, yet calls their stats helper whenever stats samples. The
  conjunction of absent optional output and requested true-correlation stats is
  undefined in Native. The draft retains source call/guard order, without
  claiming numerical equivalence for undefined input data or adding JAX's guard.
- `l_diagnose_correlations && l_calc_w_corr` reaches the lower diagnose owner's
  explicit JAX-supported exclusion: its Native working matrix is uninitialized.
  This is a lower API gate, not a blanket exclusion of `l_calc_w_corr`.
- Native sampled `rtp2_from_chi` uses `rt_2(j,:)` after the stats hydrometeor DO
  loop, so `j=hydromet_dim+1`. The draft preserves that column. Native lines
  753–780 and 4369–4375 confirm it affects the diagnostic: with two hydrometeors
  and four columns, column three is used for every column; with two hydrometeors
  and one column, that read is outside the array. JAX uses each corresponding
  column. No Native fix or replacement value is introduced. This stats branch
  has no defined BFB reference when the Native read is outside its bounds.

The first fatal negative-hydrometeor check uses an integer reduction only to
find whether failures exist. Host inspection retains Native species/level/column
order, prints the first failure, marks only that column and returns immediately.
The precipitation debug return, all-column assertion status and subsequent error
context remain in source order. Numeric diagnostic values use local
`max_digits10`; compiler-specific Fortran spacing is unverified.

All CPU zero-tolerance evidence remains pending integration, including ordinary
and nondefault precipitation/correlation options, optional output, dynamic stats,
clipping, matrix failures, threshold cases and hydrometeor counts. GPU comparison
is an occasional later milestone, not an additional draft completion gate.
