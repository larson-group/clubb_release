# JAX Conversion Plan for `src/CLUBB_core`

## Overview

This document records the architecture and maintenance workflow for the completed
JAX conversion of the supported [`src/CLUBB_core/`](./src/CLUBB_core/) path.
The active driver is computationally independent of `clubb_python_api`; the
Fortran standalone is an external regression oracle, not a runtime dependency.

Numerical equivalence with the current Fortran path remains the primary
engineering requirement. Structural mirror quality, JIT compatibility, and
explicit fail-loud boundaries for unsupported model features are also acceptance
criteria.

## Current Architecture

The supported standalone path is fully JAX-backed:

- [`clubb_jax/`](./) provides a runnable SCM path through
  `./run_scripts/run_scm.py -jax ...`.
- [`run_jax_vs_fortran_cases.py`](../tests/run_jax_vs_fortran_cases.py)
  compares the independent Fortran and JAX standalone outputs.
- `src/advance_clubb_to_end.py` owns the Python lifecycle around the jitted
  `advance_clubb_core` entry point.
- Shared model objects and statistics state are JAX-owned pytrees. NetCDF and
  other host I/O remain outside compiled physics.
- Unsupported physics and lifecycle combinations are rejected during
  initialization instead of falling back to Fortran.

## Scope of the Port

The conversion target is the supported standalone execution surface of
[`src/CLUBB_core/`](./src/CLUBB_core/), ported faithfully to JAX, together with
the driver code required by
[`clubb_jax/src/advance_clubb_to_end.py`](./src/advance_clubb_to_end.py).

This scope definition is important:

- The timestep-loop path defines the required converted surface and the primary
  regression scope.
- Inactive or unsupported Fortran branches must be documented and rejected;
  they must not be silently delegated to F2PY.
- Additional driver code may be converted where needed to support the JAX
  timestep loop or to simplify the final architecture, but driver-file parity
  with `clubb_python_driver/` is not the goal.
- Success is defined by faithful `CLUBB_core` behavior on the supported path,
  not by mechanically converting unrelated source trees.

The core conversion phase is complete. Current work should focus on regression
testing, structural audits against changed Fortran, correctness fixes, and
measured performance improvements.

## Radiation Path

The supported non-BUGS radiation path is now computationally independent of
the old mutable Python driver:

- `Radiation/radiation_module.py` owns the active functional radiation call,
  its source-order routines, and the per-sample radiation stats cadence.
- `simple_rad_module.py`, `rad_lwsw_module.py`, `cos_solar_zen_module.py`, and
  `soil_vegetation.py` are directly JIT-compatible JAX ports. LBA table I/O
  remains a host initialization boundary.
- `RadiationParameters` is the sole documented adaptation for Fortran's
  `parameters_radiation` module globals. It is an immutable JAX pytree with
  lookup tables as leaves and configuration as static metadata.
- BUGSrad and SILHS radiation remain explicit initialization-time unsupported
  boundaries. Their modules are not runtime fallbacks.

The outer timestep loop calls radiation every timestep, passing `l_rad_itime`
to preserve the source distinction between recomputing fluxes and updating
sampled diagnostics.

## Recommended Maintenance Workflow

Fortran changes should be synchronized incrementally rather than by rewriting
or replacing the JAX tree wholesale.

Use these documents in this order:

- `LLM_prompts/SHORTCUTS.md` is only the discovery index.
- This file records the whole-core architecture and acceptance criteria.
- `LLM_prompts/port_underlying_fortran_to_other_languages.md` is the executable
  structural checklist for changed files.
- File-specific plans, such as `advance_xm_wpxp_pytree_migration_plan.md`,
  are historical checkpoints and examples, not the default starting point for a
  new file.

Recommended working rules:

- Start with the canonical changed Fortran routine and audit its JAX mirror
  block-by-block.
- Update one file, routine family, or tightly coupled batch at a time.
- Keep the JAX driver runnable after each incremental step.
- Never add a `clubb_python_api` or raw F2PY fallback. Reject an unsupported
  branch explicitly until it is ported.
- Preserve behavior first. Avoid opportunistic redesigns that make correctness
  regressions harder to localize.

The exact file-by-file order should be deferred to the implementation team.
However, work selection should follow these principles:

- Prioritize routines on the active timestep path.
- Prioritize changes that affect the active timestep loop or curated cases.
- Group tightly coupled modules only when splitting them would create excessive
  interface churn.
- Avoid early restructuring for JAX style, performance, or aesthetics unless it
  clearly reduces transitional complexity.

## Type and Interface Considerations

Fortran derived types are represented by JAX-owned immutable objects. Their
field names and behavior should stay close to the source types, while array
leaves and static metadata must remain explicit for pytree/JIT behavior.

The main mirrored model-state types are:

- `grid`
- `sclr_idx_type`
- `clubb_config_flags_type`
- `nu_vertical_res_dep`
- `pdf_parameter`
- `implicit_coefs_terms`
- `err_info_type`
- `stats_type`

Object rules:

- Put repeated domain operations on the shared object interface rather than
  duplicating tiny helpers inside each ported file. For example, `ErrInfo`
  should expose functional methods such as `is_fatal()`, `set_fatal(mask=...)`,
  and `reset_code()`.
- These methods must return updated objects rather than mutating in place. This
  keeps the port surface close to the Fortran intent while preserving explicit
  state threading for JAX.
- JAX-owned state containers should be JAX-compatible pytrees, preferably
  frozen dataclass-like objects with clear update methods such as
  `state = state.replace(...)` or domain-specific equivalents.
- Host conversion belongs only at explicit I/O boundaries. Ported
  computational routines must not contain converter plumbing.

## Stats Strategy

The earlier deferred-stats phase is complete. The active JAX driver no longer
uses the F2PY stats API.

- `src/CLUBB_core/stats_netcdf.py::StatsWriter` owns registry expansion, cadence and time
  windows, optional batching/remapping/SILHS output, and the
  NetCDF schema mirrored from
  [`stats_netcdf.F90`](../src/CLUBB_core/stats_netcdf.F90).
- `src/CLUBB_core/jax_stats.py::JaxStats` is an immutable pytree whose public
  `update`, `begin_budget`, `update_budget`, `finalize_budget`,
  `var_on_stats_list`, and `l_sample` surface preserves existing core call
  sites. Diagnostic payloads use `stop_gradient`.
- Accumulators are grouped into fixed-shape banks by source grid. Compiled code
  updates those banks directly; there is no event log or per-timestep replay.
- The host receives all banks in one `jax.device_get` only at an output
  boundary, then averages and writes the record through `StatsWriter`.
- The root `run_scm.py -jax` path launches only the repository JAX virtualenv
  and does not select or add a compiled Fortran Python runtime.

The active standalone driver still rejects SILHS, adaptive-grid, radiation-grid,
and multi-batch configurations at its existing feature gates. The backend
implements those stats interfaces for future callers, but enabling the related
model subsystems remains separate driver work.

## Verification and Acceptance Criteria

The recurring regression workflow should be:

- Run [`run_jax_vs_fortran_cases.py`](../tests/run_jax_vs_fortran_cases.py)
  after each meaningful incremental conversion.
- Use `./run_scripts/run_scm.py -jax <case>` for targeted case-level debugging.
- Treat any new mismatch as a likely bug first, even if only one or two cases
  fail.
- If many cases fail immediately after one porting step, treat that as strong
  evidence that the new change introduced a defect.
- Only consider accumulated roundoff after targeted bug review does not reveal
  a concrete issue.

The comparison harness is worth trusting as a primary signal because it has
already shown robust agreement across compiler/compiler-setting changes,
Python-driver execution, and CPU/GPU contexts. Because of that, any mismatch
should be treated as actionable unless proved otherwise.

Recommended acceptance criteria:

- **Change-level success:** the JAX driver remains runnable and the comparison
  harness does not introduce unexplained new mismatches after each step.
- **Core-path success:** code exercised from `src/advance_clubb_to_end.py` has no
  `clubb_python_api` dependency and unsupported branches fail clearly.
- **Release success:** no API usage remains in the JAX timestep path, and the
  curated JAX-versus-Fortran cases pass without a compiled Python/Fortran bridge.
- **Structural success:** every changed radiation module must pass a final
  block-by-block review against
  `LLM_prompts/port_underlying_fortran_to_other_languages.md`. Numerical
  agreement does not excuse changed routine inventory, ordering, signatures,
  comments, call order, calculation order, or undocumented JAX adaptations.

## Risks and Interpretation Guidance

The main technical risks are:

- subtle interface mismatches between JAX objects and their Fortran semantics
- incorrect assumptions about pytree static/dynamic fields or array layout
- regressions that appear numerical but are actually logic bugs introduced
  during a single incremental port
- compilation regressions caused by changing static metadata or JIT boundaries

The main process risks are:

- expanding scope from “code reachable from the timestep loop” to “everything
  in the cloned driver”
- broadening a focused source synchronization into unrelated driver work
- treating small output divergence as acceptable too early, before concrete bug
  investigation has been exhausted

The recommended working posture is conservative:

- keep the executable path alive
- change one meaningful piece at a time
- validate after each step
- assume mismatches are bugs until shown otherwise
- measure JAX-specific optimizations and retain only those that preserve the
  comparison baseline
