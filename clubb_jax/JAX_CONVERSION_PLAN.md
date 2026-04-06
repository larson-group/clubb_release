# JAX Conversion Plan for `src/CLUBB_core`

## Overview

This document describes the recommended workflow for converting
[`src/CLUBB_core/`](../src/CLUBB_core/) from Fortran to JAX using the current
[`clubb_jax/`](./) driver together with
[`clubb_python_api/`](../clubb_python_api/) as a transitional scaffold.

The current `clubb_jax/` path is a working staging environment, not yet a true
JAX implementation. Its purpose is to make a top-down, incremental port
practical while preserving a runnable SCM path throughout the conversion.

The primary engineering goal during this work should be numerical equivalence to
the current Fortran path. JAX-native cleanup, JIT strategy, autodiff support,
and accelerator-oriented optimization should be deferred until correctness is
stable.

## Current Transitional Architecture

The present setup exists to support incremental replacement rather than to serve
as the final design:

- [`clubb_jax/`](./) provides a runnable SCM path through
  `./run_scripts/run_scm.py -jax ...`.
- [`clubb_python_api/`](../clubb_python_api/) exposes the public cross-module
  routine surface of `src/CLUBB_core`, which allows missing JAX pieces to stay
  Fortran-backed while neighboring logic is ported.
- [`run_jax_vs_fortran_cases.py`](../run_scripts/run_jax_vs_fortran_cases.py)
  provides an incremental comparison harness against the Fortran standalone.
- This combination allows the driver to remain operational while individual
  core files or routine families are replaced.

The intended end state is a JAX path with no
`clubb_python_api` usage inside the timestep loop, and ideally no API
dependency at all.

## Scope of the Port

The ultimate conversion target is the entirety of
[`src/CLUBB_core/`](../src/CLUBB_core/), ported faithfully to JAX, together
with any additional driver code required to run the timestep loop through
[`clubb_jax/advance_clubb_to_end.py`](./advance_clubb_to_end.py).

This scope definition is important:

- The timestep-loop path defines the critical path for incremental progress,
  because it is the portion that must become JAX-backed first in order to
  produce a computationally independent runnable driver.
- Full `CLUBB_core` conversion remains the intended end state, even for files
  not exercised immediately by the active timestep path.
- Additional driver code may be converted where needed to support the JAX
  timestep loop or to simplify the final architecture, but driver-file parity
  with `clubb_python_driver/` is not the goal.
- Success should be defined by faithful `CLUBB_core` conversion and removal of
  transitional API use from the computational timestep path, not by
  mechanically converting every copied file under `clubb_jax/`.

In other words, full `CLUBB_core` conversion is the goal, while the
`advance_clubb_to_end.py` timestep path determines the recommended order of
operations.

## Recommended Incremental Workflow

The port should proceed top-down and incrementally rather than as a full
rewrite.

Recommended working rules:

- Port from routines already exercised by
  [`advance_clubb_to_end.py`](./advance_clubb_to_end.py).
- Replace one file, routine family, or tightly coupled batch at a time.
- Keep the JAX driver runnable after each incremental step.
- Where JAX code is not ready yet, continue using `clubb_python_api`
  internally to fill that gap.
- Prefer explicit keyword arguments when calling API wrappers during the
  transition so interface mismatches are easier to audit.
- Preserve behavior first. Avoid opportunistic redesigns that make correctness
  regressions harder to localize.

The exact file-by-file order should be deferred to the implementation team.
However, work selection should follow these principles:

- Prioritize routines on the active timestep path.
- Prioritize replacements that remove API calls from inside the timestep loop.
- Group tightly coupled modules only when splitting them would create excessive
  interface churn.
- Avoid early restructuring for JAX style, performance, or aesthetics unless it
  clearly reduces transitional complexity.

[`clubb_jax/advance_clubb_core.py`](./advance_clubb_core.py) should be treated
as experimental reference material rather than the target architecture. It may
still be useful when interpreting the intent of existing Python code, but it
should not be treated as the authority for final JAX design decisions.

## Type and Interface Considerations

The transition risk around Fortran-derived-types and Python-object mirrors
should be treated explicitly.

The current API converters expect specific Python-side storage layouts for the
mirrored Fortran types. JAX-side replacements for those Python objects are
desirable, but while the API remains in use they must preserve
converter-compatible storage.

This is a compatibility constraint during the transition, not a design
preference.

The main bridge-relevant types are:

- `grid`
- `sclr_idx_type`
- `clubb_config_flags_type`
- `nu_vertical_res_dep`
- `pdf_parameter`
- `implicit_coefs_terms`
- `err_info_type`
- `stats_type`

Additional guidance:

- The API is well tested and generally reliable, so it is a good transitional
  backend.
- If results diverge, an API issue is still possible even if it is less likely
  than a JAX-port bug.
- During debugging, newly ported JAX code should be suspected first, but the
  API should not be assumed perfect.
- If JAX-native versions of the current Python-object mirrors are introduced
  before the API is removed, they should continue to respect the storage model
  expected by the existing converters.

## Stats Strategy

[`stats_netcdf.F90`](../src/CLUBB_core/stats_netcdf.F90) and `stats_type` should be
treated as a special case.

Relevant characteristics:

- `stats_type` is the only major stored type without a Python equivalent today.
- It is high complexity and high data volume.
- Moving it back and forth between Python and Fortran during runtime is costly.
- It is low priority relative to the computational timestep path because it
  primarily handles output rather than core model evolution.

Recommended default:

- Defer stats conversion until the computational port is stable and validated.

The implementation team should still be aware of the alternatives:

- Port stats alongside the core. This remains possible, but it would likely
  slow the overall effort and force difficult `stats_type` work much earlier
  than is otherwise necessary.
- Leave stats permanently Fortran-backed as an accepted backend dependency for
  output only.

The default recommendation is to avoid letting stats become the blocker for the
core conversion effort. If a later decision is made to port stats fully, that
work should follow confirmation that the computational path is already correct.

## Verification and Acceptance Criteria

The recurring regression workflow should be:

- Run [`run_jax_vs_fortran_cases.py`](../run_scripts/run_jax_vs_fortran_cases.py)
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

- **Incremental success:** the JAX driver remains runnable and the comparison
  harness does not introduce unexplained new mismatches after each step.
- **Core-path success:** code exercised from `advance_clubb_to_end.py` no longer
  depends on `clubb_python_api` for computation.
- **Final success:** no API usage remains in the JAX timestep path, and any
  remaining Fortran dependency is either eliminated or explicitly justified,
  such as a deliberate choice to keep stats/output Fortran-backed.

## Risks and Interpretation Guidance

The main technical risks are:

- subtle interface mismatches between JAX code and the existing API boundary
- incorrect assumptions about mirrored-type storage layout
- regressions that appear numerical but are actually logic bugs introduced
  during a single incremental port
- premature optimization or JAX-oriented redesign before the behavior is
  stable

The main process risks are:

- expanding scope from “code reachable from the timestep loop” to “everything
  in the cloned driver”
- allowing stats/output complexity to dominate the early schedule
- treating small output divergence as acceptable too early, before concrete bug
  investigation has been exhausted

The recommended working posture is conservative:

- keep the executable path alive
- change one meaningful piece at a time
- validate after each step
- assume mismatches are bugs until shown otherwise
- defer JAX-specific optimization goals until correctness is secure
