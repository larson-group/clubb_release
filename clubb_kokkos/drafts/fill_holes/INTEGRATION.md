# Fill negative fields while conserving available mass

This isolated draft covers all eight Native public `fill_holes.F90` routines,
plus its private global and sliding-window methods. Complete Native and actual
JAX scope were read before the source outline and implementation. No production
source or shared generator was edited. Compilation and numerical checks are pending.

Use the Native-named header/source as one module when root reaches this live
boundary. Public argument names/order, InOut/Out ownership, optional five stats
names and Native one-based grid/species indices are retained. Batched physical
fields use column-first LayoutLeft. Single-column microphysics interfaces retain
their Native `(level,species)` axes; the driver packs strided column sections
through the existing `array_types.hpp` owner and restores caller allocations.

Global filling preserves each column's original numerator, denominator, clipped
numerator and eps/2 guard. Sliding windows retain source-direction overlapping
updates, denominator sums, all-field hole checks and global fallback. These
integrals/recurrences run sequentially within independent column kernels; no
parallel floating-point reduction is substituted. Same-phase redistribution
packs only the original frozen mixing-ratio subset in species order. Its
device-callable three-argument overload lets all levels execute in one kernel;
the concrete public wrapper and batched parent render diagnostics on the host.

`fill_holes_vertical_api` supports global and sliding-window methods, as actual
JAX does (`fill_holes.py:1–13,235–323`). Native widening windows, smart window,
smooth smart window and parallel fill are explicitly outside this supported
scope; no fallback silently substitutes for them. Dispatch still happens only
after finding a hole. Both grid directions use Native bounds directly rather
than JAX's axis reversal. No new dimension or division guards were added.

Two public lower metadata calls are pending in the existing f2k owner:

| Native `index_mapping` function | Required input | Return |
| --- | --- | --- |
| `Nx2rx_hm_idx` | `Nx_idx` integer, `hm_metadata` In UDT | integer `rx_idx`, Native one based |
| `mvr_hm_max` | `hydromet_idx` integer, `hm_metadata` In UDT | real `mvr_hydromet_max` [m] |

Both execute on the host once per species, outside all device kernels. The
declaration-only `f2k_dependencies.hpp` records these actual public contracts;
replace it with generated definitions at promotion. No raw binding, global
Fortran storage or scalar device crossing is introduced. Add `global_fill=1`
and `sliding_window=2` through the existing generated `model_flags` owner.
Other required constants and stats begin/finalize APIs already exist.

Native budgets preserve name/species/stage order and sample gates, including
blank names. A local driver lambda stages `field/dt` immediately at the host
stats call. Host warning passes retain old/new values, source order and partial
fill warning multiplicity. Unqualified invalid-selector `STOP` maps to normal
process exit; unknown vapor species reaches error exit only upon an actual hole
and after the Native vapor-tendency mutation. Text values use C++ rendering;
Fortran list-directed formatting is not claimed byte-identical.

Source review includes the original unfinished liquid same-phase filling and
unsupported concentration mappings: no guessed species or radius fallback is
added. GPU performance of per-column microphysics packing is unmeasured; this
is an explicit public API baseline, not a fundamental GPU requirement.

Run the maintained CPU whole-model zero-tolerance comparison after hookup.
Descending grids, optional stats names, debug warning branches, invalid method
termination and unsupported vapor-species exits currently have source review
only. GPU execution is a later checkpoint under the current validation policy.

No deliberate FMA, forced-rounding or math workaround was introduced for BFB.
Integer powers use direct multiplication, EPSILON uses the real kind, and
source association/order remain. Numerical preservation is unverified.
