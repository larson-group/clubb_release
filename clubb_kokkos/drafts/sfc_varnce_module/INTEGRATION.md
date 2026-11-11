# Surface-variance draft

[sfc_varnce_module.cpp](sfc_varnce_module.cpp) ports the complete 33-argument
`calc_sfc_varnce` procedure from
[Native sfc_varnce_module.F90](../../../src/CLUBB_core/sfc_varnce_module.F90).
The [header](sfc_varnce_module.hpp) preserves the Native argument order and ranks.
The full source was read and [outlined](SOURCE_OUTLINE.json) before code.

**Compilation, integration, execution and numerical/BFB validation are pending.**
This isolated folder changes no production source, CMake input or generator.

## Promote with its reached caller

Copy the `.cpp`/`.hpp` into `clubb_kokkos/src/CLUBB_core/`, register the source
through the existing CMake owner, and replace the reached whole-procedure
`f2k::calc_sfc_varnce` call. Remove that obsolete bridge contract at promotion.
No new physical/model constant is required.

The existing Native `stats_begin_budget` and `stats_finalize_budget` interfaces
remain in use. Each sampled call materializes the Native `field / dt` array
expression into a local execution-space View immediately before that call. All
six begin calls and six finalize calls retain Native field and call order:
`thlp2`, `rtp2`, `rtpthlp`, `up2`, `vp2`, `wp2`. These buffers are C++ array-expression
storage, independent of the chosen stats backend. They introduce no persistent
bridge state or new stats API.

One lower public diagnostic dependency is explicit:

```text
numerical_check::sfc_varnce_check(
  sclr_dim,
  wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc,
  sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc,
  err_info)
```

All 11 arguments are required: one integer, six core-real scalars, three
core-real rank-one arrays of `sclr_dim`, and mutable `err_info_type`.
[expression_audit.json](expression_audit.json) records the complete parsed Native
contract. The draft calls its public `f2k::sfc_varnce_check` name, which the current
production bridge has not yet declared. Add that temporary public contract via
the existing generator owner, or promote the corresponding numerical-check
routine and replace this call directly. No whole surface procedure is delegated.

Do not substitute JAX's inline checker: Native checks non-finite values,
including infinities, and a failed supplied value can mark the entire
`err_info.err_code` array fatal. The draft stages the nine checked fields once
inside the Native debug-level-2 gate, calls the checker by column, prints headers
at the same point and retains the final complete-array input/output dump order.
No Fortran call occurs inside a device lambda.

## Preserved branches and boundaries

Both `l_vary_convect_depth` choices and passive-scalar emulation branches are
ported. Scalar/grid metadata retain Native one-based values; View accesses
subtract one. Independent surface columns and column/scalar pairs run in
parallel. Local scratch belongs to the call, and device bodies capture only
component Views or scalar metadata.

JAX omits the source-constant-false Andre et al. (1978) branch. This complete
Native draft retains that entire branch and the `false` parameter, without
adding a runtime switch. Its wind rotation, stability powers, scalar covariance
and splatting formulas stay in source order. The active scheme preserves
variance floors, splatting corrections, `wp2_max` clipping and the later reset
when the lowest grid level lies away from the surface.

The positive-heat-flux depth search advances with `grid_dir_indx` until the flux
is nonpositive or height reaches 1000 m. Native does not include an upper-grid
bound; the draft preserves those stopping criteria and does not import JAX's
last-level fallback. Input profiles therefore retain the Native stopping
assumption. Fractional powers and other dormant-branch domain assumptions are
also unchanged.

MIN/MAX selectors retain the established second-tie convention (`BFB-003` in
[the adaptation log](../../BFB_ADAPTATIONS.md)); three-argument MAX is expanded
in argument order. Integer powers use ordinary grouped products and real powers
use `Kokkos::pow`. No explicit FMA, special-function emulation or speculative
rounding workaround was added.

## Review and runtime evidence

[SOURCE_AUDIT.json](SOURCE_AUDIT.json) records the separate source cleanup and
manual review. [source_audit.py](source_audit.py) compares the 33-argument
contract, all 80 science assignment trees, 24 ordered branch predicates, all
12 sampled budget field/call mappings and ownership of 36 device bodies. It is
a source review aid; it does not compile or execute physics. The checker staging,
column/error/dump order, while-loop bounds and comments were reviewed separately.

Run the maintained zero-tolerance CPU Kokkos/CPU Fortran model comparison after
promotion. GPU comparisons are occasional milestones; record their differences
without making GPU BFB a per-port prerequisite. No model/test run is claimed here.
