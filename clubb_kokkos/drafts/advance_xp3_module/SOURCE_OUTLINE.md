# Third-moment advance source outline

Read all 1006 native lines and the complete actual JAX module before porting.
Native source is the authority. JAX declares only state/stats representation
adaptations; it excludes no native scientific branch in this module.

| Native routine, source order | Lines | Scope |
| --- | --- | --- |
| advance_xp3 | 36–283 | Public predictive/steady-state third-moment parent |
| diagnose_xp3 | 286–566 | Public diagnostic third-moment parent |
| advance_xp3_simplified | 569–856 | Private scalar/field update and budget terms |
| term_tp_rhs | 859–933 | Private device scalar production expression |
| term_ac_rhs | 936–1002 | Private device scalar accumulation expression |

All native arguments are required, with no optionals. Exact inventories are in
native_inventory.json. Public inputs/outputs follow native order and units.

advance_xp3: interpolate wp2 with its minimum; call simplified updates in order
rt, thl, passive scalars; compute vertical skewness; interpolate scalar/wind
fluxes and variances with their original distinct minima; compute moist static
stability and xp3 coefficient; diagnose u/v third moments; sample eight original
interpolated fields. Private scalar calls keep their source loop/call order.

diagnose_xp3: same skewness and field interpolation sequence; branch on ADG1.
ADG1 uses interpolated PDF width for thl/rt/u/v/scalar ansatz calls. Other PDFs
compute static stability, explicitly zero the unused width field, compute the
modified coefficient, then call the same thl/rt/u/v/scalar ansatz in native order.
Preserve the ADG1 passive scalar flux's native sclr_tol**2 interpolation minimum;
the non-ADG1 flux interpolation has no such minimum. Sample the same eight fields.

advance_xp3_simplified: select rt/thl/empty scalar budget names; retain the local
l_predict_xp3=false and C_xp3_dissipation=1 parameters with both predictive and
steady-state branches. Initialize production/accumulation arrays; interpolate
mean and variance; independent column/level cells compute the two scalar helpers,
then predictive/LMM or steady-state update. Preserve source k=1:nzt-1 and
kp1=min(k+1,nzt), even though the nearby comment says km1. Zero storage level nzt
as the native upper boundary for both grid directions. Sample production,
accumulation and dissipation; retain predictive-only begin/finalize budgets.

The two scalar helpers preserve the complete native scientific discretization
comments and operation associations. They must be device callable inside the
cell loop. No vertical recurrence or matrix solve exists in this module.

Existing live APIs cover all lower requirements: grid zm2zt_api/zt2zm_api/ddzm,
Skx_func/xp3_LG_2005_ansatz, f2k stats update/begin/finalize, core constants and
parameter/PDF indices. No missing lower declarations or new bridge is planned.
Use existing array_types staging only where a scalar/subview must cross those
typed boundaries. No host I/O, fatal/error path, random draw or mutable global
state exists in the native file.

Implement in an isolated source/header, then perform separate signature/call/
formula/branch/comment/format audit. No tests or builds. Root owns live hookup,
regular CPU/CPU zero-tolerance model validation and occasional GPU checkpoints.
No speculative FMA or other bit workaround is authorized by this draft task.
