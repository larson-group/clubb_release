# Cloudy updraft diagnostic draft

`pdf_closure_module.cpp/.hpp` contains the next public Native routine,
`calc_w_up_in_cloud`, with all 13 required arguments. It computes four
caller-owned output fields in one column-first cell kernel. It has no lower
physics calls, bridge declarations, allocations, global storage, statistics or
error side effects. The existing constants owner supplies every required
constant.

Merge the definition after the private `calc_xprcp_component` and before
`pdf_closure_driver`, outside the unnamed namespace, preserving Native order.
Add its declaration in the same public position and replace the current
`f2k::calc_w_up_in_cloud` call in `pdf_closure`. Remove that obsolete bridge
through its existing owner. These files are fragments for the permanent
`pdf_closure_module` owner; they are not a separate production module or target.

Keep the caller's current ADG1/ADG2/new-hybrid and requested-statistics gates.
Native calls this routine only when either velocity diagnostic is requested;
the existing caller writes all four fields to zero otherwise. No new PDF
argument or inner option guard belongs in this helper.

The implementation retains both components' three distinct branches, Native
scalar names, strict thresholds, mixed downdraft fraction `one-updraft_frac`,
and the four final output expressions in source order. Integer squares use
ordinary grouped multiplication under the original unary minus. Conditional
MAX lowerings retain the existing second-operand tie choice, including its
unordered-comparison behavior; there is no `fmax` replacement or extra floor.
There is no explicit FMA, alternate math library or measured BFB adjustment.

The two `varnce_w_*` comments say **variance** in the target. Native calls them
standard deviations despite `[m^2/s^2]` units and immediately takes their square
roots to obtain `stdev_w_*`. This corrects wording only.

`SOURCE_OUTLINE.json` was written before implementation. `audit_source.py` and
`source_audit.json` record a separate signature/branch/operation-tree/comment
review. Cleanup checked source groups, local ownership, all branch writes and
formatting. No compilation, test or numerical execution was performed. Root
owns integration and the maintained CPU Fortran/CPU Kokkos zero-tolerance gate;
GPU comparisons remain occasional milestones. Source agreement does not
establish bitwise equality of `sqrt`, `erf`, `exp` or the full output formulas.

## Unreachable interpolation cluster

`unused_interpolation_scope.json` records all 14 symbol occurrences found by a
case-insensitive search of **all `src/` and `clubb_jax/`**. Every occurrence lies
inside the three-function Native cluster at lines 3343–3589:

- `interp_var_array` calls `var_value_integer_height` and `var_subgrid_interp`.
- The latter two call the public lower `mono_cubic_interp`.
- No source outside the cluster references any of the three names; no entry
  calls `interp_var_array`. The Native module defaults to private and exports
  only its seven listed public procedures, excluding all three helpers.
- The three module import sites under `src/` select only `pdf_closure_driver`,
  `calc_wp4_pdf` or `calc_w_up_in_cloud`. No generic interface or procedure
  pointer introduces another named edge.
- Actual JAX defines and references none of the three helpers anywhere in its
  tree. Its cloudy-updraft helper contains the live computation instead.

Excluding this private, source-unreachable cluster therefore matches the
current JAX port's actual scope. This is a static current-repository conclusion;
it neither removes Native source nor claims runtime coverage. Revisit the
exclusion if a future Native entry point adds an incoming call.
