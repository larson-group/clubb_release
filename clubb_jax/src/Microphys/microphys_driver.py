"""JAX port of microphys_driver.F90 — the per-step microphysics scheme dispatch.

Mirrors clubb_release/src/Microphys/microphys_driver.F90: `calc_microphys_scheme_tendcies` — dispatch the
configured per-timestep microphysics scheme (Khairoutdinov-Kogan or Morrison M2005), computing + applying its
tendencies to the prognostic state. The per-scheme implementations live in their own step modules
(Microphys/kk_microphys_step.py, Microphys/morrison_microphys_step.py); advance_clubb_to_end's timestep loop calls
this dispatch, mirroring the Fortran driver's `call calc_microphys_scheme_tendcies`.
"""

from __future__ import annotations


def calc_microphys_scheme_tendcies(state: dict, time_current: float) -> None:
    """Dispatch the per-step microphysics scheme (microphys_driver.F90:calc_microphys_scheme_tendcies).

    A no-op for `microphys_scheme='none'`. For Morrison it is also skipped until `microphys_start_time`
    (microphys_driver.F90:389 — e.g. nov11 has a 60-step spin-up; before it, no `*_mc` and no hydrometeor advance).
    """
    scheme = state.get('microphys_scheme', 'none')

    # ── KK (Khairoutdinov-Kogan) rain microphysics ──────────────────────────
    # Computes + stores the KK tendencies on live state, then (gated behind l_kk_micro_apply, which clubb_driver
    # sets TRUE for KK cases — clubb_driver.py:1069) applies the hydrometeor transport + feedback. So the full
    # KK path is live for KK cases (rico, BLOCKED/FP-limited); see kk_microphys_step.py and DESIGN.md "rico".
    if scheme == 'khairoutdinov_kogan':
        from clubb_jax.src.Microphys.kk_microphys_step import advance_kk_microphysics
        advance_kk_microphysics(state)

    # ── Morrison (M2005 2-moment) microphysics ──────────────────────────────
    # Computes the CLUBB-form *_mc (rcm/rvm/thlm + 8 hydrometeors) via morrison_microphys_driver and advances
    # the hydrometeor fields.
    elif scheme == 'morrison' and time_current >= state.get('microphys_start_time', 0.0):
        from clubb_jax.src.Microphys.morrison_microphys_step import advance_morrison_microphysics
        advance_morrison_microphysics(state)
