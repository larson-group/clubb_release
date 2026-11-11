# Scalar-variance timestep draft

The complete supported `advance_xp2_xpyp_module.F90` bodies are translated here:
all 16 Native routines, including the 66-argument public timestep, 18-argument
rain-evaporation update and fourteen private procedures. The actual single- and
multiple-LHS solvers are present; JAX's unfinished top-level rejection of the
multiple-LHS route is not copied.

**Compilation, integration, runtime, numerical/BFB and performance evidence are
pending.** No production, CMake, shared generator, tests or root worklog were
changed. Root owns promotion and CPU Fortran versus CPU Kokkos BFB at zero
tolerance before integration advances. GPU runs provide occasional milestone
evidence, with differences documented.

The material scope exclusion is the Native up2/vp2 sponge block, matching JAX's
upstream rejection of sponge-enabled cases. The exact source boundary is marked
in the C++ file. Promotion must enforce that rejection before this routine; it
must not accept a sponge case and silently omit its damping. All other source
routines, field-specific budgets, clipping, LMM, Coriolis, scalar emulation,
centered/upwind/Godunov branches and fatal-return paths remain present.

Completed separate source audit and cleanup:

- Read all 6341 Native lines and wrote the full routine/argument/branch outline
  before implementation; inspected JAX for supported scope and exclusions.
- Verified all 16 Native signatures and relative definition order, both public
  headers, and 1864 non-solver call-argument mappings in source order.
- Independently compared 345 explicit assignment operation trees, 85 whole-array
  assignment trees with their slice indices, 175 ordered IF predicates and all
  non-solver CASE identifier orders. MIN/MAX and integer-power grouping are
  included; no math execution is used by these checks.
- Checked interpolation function call arguments separately, all original
  diagnostic write labels/order, 257 capture-safe device bodies and the
  documented matrix-layout, reduction and ordered-recurrence adaptations.
- Re-read the solver's condition-number choices, raw InOut restoration, Native
  fatal return before final solution reversal, clipping/statistics field order,
  and local storage ownership. Applied two-space/100-column C++ formatting and
  checked delimiters/trailing whitespace.

[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json) records the Native contracts and initial
scope decisions. [source_audit.json](source_audit.json) and
[formula_audit.json](formula_audit.json) record reviewed source hashes and the
source-only checks. The audit includes declarations of legacy missing INTENT
contracts verified from their complete bodies; Native was not modified.

Run `python3 clubb_kokkos/drafts/advance_xp2_xpyp_module/audit_formulas.py` to
recreate both source-review records. This is a source audit, not a numerical
test. `prepare.py` recreates the initial draft and must be followed by formatting
and audit; reconcile later numerical fixes before rerunning it.
[INTEGRATION.md](INTEGRATION.md) lists the unresolved public lower contracts and
[BFB_NOTES.md](BFB_NOTES.md) separates ordinary translation from pending evidence.
