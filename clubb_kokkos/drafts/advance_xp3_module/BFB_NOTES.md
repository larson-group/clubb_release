# Third-moment arithmetic notes

No forced FMA, reassociation, alternate math approximation, new clipping guard
or other speculative bit workaround was introduced. This draft has no compiled
or numerical evidence. CPU Fortran versus CPU Kokkos at zero tolerance is the
regular integration gate; GPU comparisons are occasional checkpoints.

Ordinary translations preserve native source operation trees:

- Source tolerance squares use operand multiplication. The xp3 coefficient
  retains its source addition/product/exp grouping and has no new clipping.
- Buoyancy-frequency MAX retains its two operand order; an ordinary local
  `stability` evaluates the first operand once. The steady-state MIN similarly
  evaluates its reciprocal once and selects the second operand on a tie.
  These are not measured cross-compiler claims for NaNs.
- The scalar production expression remains
  `+three*xp2_zt*invrs_rho_ds_zt*invrs_dzt*(rho_ds_zmp1*wpxpp1-rho_ds_zm*wpxp)`.
  Accumulation remains `-three*wpxp2*invrs_dzt*(xm_zmp1-xm_zm)`.
- The predictive update retains the native sum/division grouping; LMM averaging
  follows the update. Steady-state multiplication/division order is unchanged.
  Dissipation statistics retain `-(C_xp3_dissipation*invrs_tau_zt)*xp3`.

Inherited source choices retained for review:

- ADG1 scalar-flux interpolation has a `sclr_tol(sclr)**2` minimum, just like
  scalar variance. Non-ADG1 scalar flux interpolation omits that minimum.
- The upper boundary is storage level `nzt` for both grid directions. The
  interior ends at `nzt-1`, and the exact native neighbor clamp is retained.
- `l_predict_xp3=false` is a private native parameter. Predictive/LMM code is
  present, but its branch has no runtime evidence from this draft.

There is no new source-definedness deviation in this module. Scratch ownership,
scalar-section adapters, integer index conversion and conditional allocation
are execution adaptations documented in INTEGRATION.md. Matching nine expression
trees and 425 call arguments establishes source structure, not BFB results.
Root must log any later deliberate arithmetic fix with exact expression,
observed cause and before/after evidence in `clubb_kokkos/BFB_ADAPTATIONS.md`.
