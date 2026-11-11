# advance_xm_wpxp_module source outline

Read the full 6472-line native source and actual JAX scope before implementation.
Only advance_xm_wpxp is public; the following 17 procedures remain private.

| Routine | Native line | Arguments |
| --- | --- | --- |
| `advance_xm_wpxp` | 59 | 107 |
| `xm_wpxp_lhs` | 1174 | 17 |
| `calc_xm_wpxp_lhs_terms` | 1429 | 27 |
| `xm_wpxp_rhs` | 1638 | 19 |
| `calc_xm_wpxp_ta_terms` | 1991 | 30 |
| `solve_xm_wpxp_with_single_lhs` | 2560 | 98 |
| `solve_xm_wpxp_with_multiple_lhs` | 3697 | 68 |
| `xm_wpxp_solve` | 4267 | 12 |
| `xm_wpxp_clipping_and_stats` | 4391 | 44 |
| `xm_term_ta_lhs` | 5190 | 7 |
| `wpxp_term_tp_lhs` | 5295 | 5 |
| `wpxp_terms_ac_pr2_lhs` | 5408 | 8 |
| `wpxp_term_pr1_lhs` | 5527 | 11 |
| `wpxp_terms_bp_pr3_rhs` | 5661 | 7 |
| `xm_correction_wpxp_cl` | 5745 | 9 |
| `damp_coefficient` | 5969 | 10 |
| `diagnose_upxp` | 6031 | 12 |
| `error_prints_xm_wpxp` | 6120 | 82 |

The exact signature and call sequence inventory is in native_inventory.json.
Only xm_wpxp_solve.rcond and xm_wpxp_clipping_and_stats.wpxp_cl_num are optional.

Parent sequence: skewness and smoothed wp3/wp2; limiter flag validation; RHS count;
LMM snapshots; C6/C7 skewness and altitude damping; coefficient stats/validation;
eddy diffusion; monotonic range; pressure coefficient; PDF advection terms; common
LHS terms; unique/equal LHS solve branch; LMM average; full fatal field dump;
excluded sponge boundaries; prognostic-wind nudging and reference stats.

Private sequence: explicit native interleaved means/flux matrix and RHS assembly;
PDF explicit/ADG1/new/new-hybrid advection variants; equal-LHS packed multi-RHS solve
with Coriolis, perturbation and horizontal scalar flux diagnoses; separate native
multiple-LHS solves; optional descending permutation and matrix API; implicit stats,
monotonic limiter, positive-definite/hole fill/relaxed covariance clipping and budgets;
individual production/pressure/buoyancy coefficients; clipping correction; damping;
horizontal flux diagnosis; complete source field diagnostics.

Retain all native calculations, option branches, stats names and fatal gates. The 71
supported stats calls keep their source position and use the Fortran backend;
eight of the native 79 calls belong to the excluded sponge blocks. JAX
explicitly excludes sponge damping because Python/JAX case initialization rejects it;
the C++ draft must mark these source boundaries and require the same upstream guard.
JAX diagnostic reductions do not constrain this draft: retain the full native dumps.

Execution plan: retain device arrays across adjacent calls, run independent source
loops on columns/cells, retain explicit host call/field/scalar order and ordered
per-column recurrences. Native bands-first matrices need an explicit permutation at
the existing matrix_solver_wrapper boundary, whose established C++ matrix layout is
(column,band,row). Preserve native-specific field and tolerance mappings.

Pending owner contracts: calc_wp3_on_wp2; calc_turb_adv_range; diffusion_zm_lhs and
diffusion_zt_lhs; term_ma_zm_lhs and term_ma_zt_lhs; xpyp_term_ta_pdf_lhs, its Godunov
form and RHS; monotonic_turbulent_flux_limit; pos_definite_adj; clip_covar;
fill_holes_vertical_api; stats_update_budget. Existing grid, Skx and matrix APIs
are used where complete. Do not edit production owners in this draft.

Source-definedness edges: xm_correction_wpxp_cl has uninitialized per-column mask;
JAX sets it false. Native positive-definite stats outputs are unwritten when no
adjustment occurs. xm_wpxp_solve.old_solution has a declaration shape inconsistent
with its callers; the excluded BiCGStab path is the only indexed use. Keep these
distinctions explicit and make no Native BFB claim for undefined values.

No speculative FMA, reduction reorder or transcendental substitutions are planned.
Source/signature/branch/comment/format audit is separate from pending compilation
and maintained whole-model numerical validation, owned by root.
