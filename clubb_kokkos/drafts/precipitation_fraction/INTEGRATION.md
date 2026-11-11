This draft ports all four routines in `precipitation_fraction.F90`, including
both component methods and the complete private assertion checker. The public
22-argument `precip_fraction` interface preserves Native argument order and
caller-owned outputs. The Native immutable public selector remains
`precip_frac_calc_type=2`; method 1 and its full weighted logic are retained
without adding a runtime option.

Promote this source/header alongside the setup PDF and Nc/Ncn owners. Include
`precipitation_fraction.hpp` and replace setup's `f2k::precip_fraction` call with
`precip_fraction` using the same 22 arguments. Delete the transient declaration
when integrating that replacement. The other dependencies are existing grid,
array/type, physical/tunable constant and error owners. Stats and the debug query
remain in their existing Native backend. No additional scientific bridge or
shared mutable storage is required.

Fields use column-first `LayoutLeft`. The public and three private section
adapters reuse the existing array owner, retain caller allocations and restore
all outputs. Logical species flags use the existing byte View representation.
The profile checker retains its Native `(level,species)` hydrometeor slice;
that local profile has no column axis. Its array adapter packs strided inputs
without changing these dimension semantics.

The two overall-fraction and two weighted-component maximum recurrences start
at physical `gr.k_ub_zt` and advance by `-gr.grid_dir_indx` to `gr.k_lb_zt`.
Each column owns its serial recurrence; levels are never launched independently
across a dependency. The column-local MAXVAL scans use ascending storage levels
and the Native negative-HUGE identity. Only immutable frozen-species flags use
integer reduction. Species boosts remain in Native species order. Metadata
flags and debug-only scalar arguments use call-local host snapshots.

`ANY(hydromet>=tol)` and `ALL(hydromet<tol)` are computed separately. They retain
their different empty-species and NaN truth values; the latter is not replaced
with the complement of the former. No JAX safe denominator, `fp==tol` identity
or optimization barrier replaces the parenthesized Native formulas. Compiler
behavior for MAX/MAXVAL ties or nonfinite inputs has not been validated.

The checker runs the original ascending-level, ordered checks in one device
work item. It returns only the first failure's integer kind/level and four field
values. Host code prints the source diagnostics, then marks **every** column in
`err_info.err_code` fatal and returns, matching the actual Native assignment.
The parent continues its source column loop and prints each header whenever
that column is already fatal. Case 6 deliberately omits the global header just
as Native does. The final mixture predicate is the one-sided residual `>eps`,
without ABS or an extra finite-value test. Diagnostic numeric values use local
`max_digits10`; exact Fortran spacing is unverified. The sampled rank-1 stats
call selects the existing explicit const-element overload.

The complete 1,180-line Native and 843-line JAX files were read before the
outline. JAX omits the debug checker from its jitted parent and provides a
reduced boolean checker; this C++ draft preserves the Native error object and
all side effects. `SOURCE_AUDIT.json` records the owner review of all 59
arguments, 62 ordered field assignments, 45 predicates, four calls with 40
actuals, and 26 device bodies. `CLEANUP_AUDIT.json` records the separate source
argument/comment/format pass. Independent review, compilation and every runtime
branch remain pending; no CPU/GPU or BFB claim is made.
