# Correlation diagnostics draft

This isolated header implements all ten active Native procedures in definition
order, including three private helpers and seven public routines. The four
wholly commented-out procedures (`approx_w_corr`, `approx_w_covar`, `set_w_corr`
and `unpack_correlations`) are excluded, matching actual JAX scope. The three
public scalar functions are device inline; array routines use concrete core
Views with four thin public generic adapters through the existing array owner.
No new lower scientific routine, metadata type, accessor or selected constant
is needed. Existing host debug-level queries remain explicit.

Promote the header before replacing the corresponding lower calls. Private
functions live in `diagnose_correlations_detail`, with forward declarations
separate from definitions. `iiPDF_w` remains Native one-based, unlike the JAX
metadata object; callers must pass the existing C++ metadata value directly.
Physical matrix `(row,column)` slots and first-index-fast LayoutLeft are retained.
All scratch Views are local to each call. Host/strided public callers stage once
through the shared owner and retain their original output allocations.

The supported `diagnose_correlations` path requires `l_calc_w_corr=false`.
Native lines87–88 initialize `corr_array_swapped` only under `.not.l_calc_w_corr`,
then line91 passes that working matrix to `diagnose_corr`, which reads its first
column at Native line163. The old approximate-w routines are entirely commented
out. JAX rejects true at its `diagnose_correlations` API entry. The draft has the
same explicit API guard before allocating/reading scratch, with a local TODO;
it does not fabricate an initial matrix.

The exact Native upstream gate is `setup_clubb_pdf_params.F90:614`
(`l_diagnose_correlations`), followed by four cloud/below/component calls at
621/626/633/638 that forward `l_calc_w_corr`. The live C++ driver currently reads
and forwards that option (`src/clubb_driver.cpp:270,615,706,1264`); it does not
supply a new gate for this draft. On integration, keep the API rejection at the
reached diagnose call, or reject the conjunction
`l_diagnose_correlations && l_calc_w_corr` at the upstream owner. A blanket
rejection of `l_calc_w_corr` would exclude other existing covariance paths that
do not reach this undefined diagnose branch. No driver edit was made here.

The false path preserves the zero scratch, array-expression SQRT argument,
rearrangement/copy/call order, ordinary Native square roots and two conditional
correlation clamps. JAX gradient floors, symmetrization and boolean-only assertion
interfaces are not imported. Disabled alternative formulas remain documented.
Three scalar mixture/correlation functions retain their exact Native expression
association and core-kind values, with integer squares expanded to products.

Approximate factorization keeps zero/one initialization, actual Native `(j,i)`
slots and ascending scalar products within each output entry. The matrix product
`MATMUL(L,TRANSPOSE(L))` uses an explicit ascending-k contraction per output cell,
with no parallel floating reduction or triangular pruning. This is ordinary
source translation; the Native compiler's intrinsic contraction/reduction order
has not been observed, so exact arithmetic equality is pending model integration.
No speculative FMA or forced-association workaround is added.

Rearrangement saves the first column and retains each of the four subsequent
Native array assignments as a separate ordered kernel. The first reversed
assignment overlaps only its unchanged final cell; all other reads remain
outside the stage's writes. Later assignments use the owned saved column.
The implementation neither substitutes a full row permutation nor discards
Native reversal order. Comments distinguish actual matrix slots where Native
upper/lower wording is misleading; call-site intents follow the callee's Out
declaration where Native copied InOut labels.

Assertions keep debug1 off-diagonal range checks before debug2 unit-diagonal
checks, their Native tolerances, exact ERROR STOP strings and check order.
ERROR STOP becomes the established C++ fatal exception boundary. They do not
add an `err_info` parameter, a finite-value check or a reduced boolean API.
No stats backend is involved.

`SOURCE_OUTLINE.json`, `SOURCE_AUDIT.json` and `CLEANUP_AUDIT.json` record source
and cleanup evidence. Compilation, all CPU/GPU/model outputs, source intrinsic
MATMUL arithmetic, swapped layouts, clipping/NaN/debug errors and host/strided
callers are untested. The true diagnose path remains explicitly unsupported.
Parent integration owns maintained CPU zero-tolerance comparisons and occasional
GPU milestones; source review is not a runtime or BFB claim.
