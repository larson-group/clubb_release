# ADG1 / ADG2 / Luhar draft

**Status: all nine native routine bodies implemented and source-audited; not integrated, compiled or numerically tested.** This folder is an independently prepared draft. No result here establishes model equivalence, Native BFB or GPU performance.

Authority is the complete current1537-line `src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90`. `SOURCE_OUTLINE.md` records the pre-implementation inventory/branches. `source_audit.json` records the final full inventory, argument/call checks, reviewed branches, intentional adaptations and source hashes. The full draft exports the same six public routines, retaining the three native private helpers and all source definition order.

## What is implemented

- ADG1 driver, w closure and private responder execute their cell calculations in column-first Kokkos Views, preserving native all-column tolerance failures and source fatal gates.
- ADG2 retains its column loop, Luhar parameter/closure calls, sigma overwrite, rt/thl and scalar responders.
- 3-D Luhar retains its per-level three-way setter choice, covariance/tolerance mappings and backsolve/closure call order.
- Luhar parameter/close functions retain both variance branches and full formulas.
- Backsolve retains the active cubic and inactive native quadratic choices, near-half/zero-skewness branches and Native HUGE value.
- Largest-root selection retains cubic/quadratic thresholds, real/complex selection and unguarded linear fallback. Native public complex root solvers are the only missing lower physics integration.

The three profile APIs currently have host C++ implementations plus explicit View-staging overloads. This avoids mixing pending native complex-root calls inside GPU kernels and keeps the formulas/control flow reviewable. Only the 3-D driver actually calls the complex-root backend. The pure calc_Luhar_params/close_Luhar_pdf routines could move to GPU independently; sharing the current host profile path is a draft implementation choice, not a fundamental requirement of those equations. It is a temporary performance boundary, not a tuned final GPU implementation. Profile adapters preserve output owners, source argument order and strided column/sequence-associated callers.

## Known source-definedness deviation

At native1022–1060, `close_Luhar_pdf` initializes local `sgn_wpxp` only under `xp2>x_tol_sqd`, but reads it in the nonvarying normalized means too. The draft initializes `sgn_wpxp = wpxp>=0 ? +1 : -1` before either branch, matching actual JAX `close_Luhar_pdf`. This choice was explicitly accepted by the root after identifying the source use.

Physical means/variances retain native formulas. No Native BFB claim applies to the native branch's undefined normalized outputs. This is a definedness deviation, not a verified BFB arithmetic workaround. Native source remains untouched; root owns any separate native follow-up/backlog work.

## Validation and remaining work

No build, test, benchmark or profile was run, as instructed. All added runtime branches are untested, including ADG2/Luhar/profile staging, w-closure GPU kernels, complex solver boundaries and inactive quadratic backsolve. The compact ADG1 parent/responder alone remains separately recoverable in `incremental_adg1/` and `/tmp/clubb_adg1_next/` for the root's incremental hookup.

Root must supply the two complex adapters declared in `f2k_dependencies.hpp`, then integrate/select appropriate source portions and run the maintained model comparison. Existing standalone Native gates require input fields for ADG2/Luhar, so their source coverage is distinct from the current19-case ADG1 matrix. Keep the native supported-option gates; source implementation does not authorize broader driver-option integration.

No source-science fixes, shared type/generator/CMake edits or production ports were made outside this owned draft folder. `expand.py` owns the full draft generation; the prior compact generator/audit/artifact are kept with the recoverable compact source. `audit_source.py` refreshes full-draft source evidence only.
