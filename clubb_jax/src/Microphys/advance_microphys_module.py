"""JAX port of advance_microphys_module.F90 — the hydrometeor transport solve.

So far: the implicit LHS contributions of the hydrometeor predictive equation that are NEW physics
(the mean + turbulent sedimentation terms). The turbulent-advection and mean-advection LHS reuse the
already-bit-faithful `diffusion_zt_lhs` and `term_ma_zt_lhs`; the sedimentation terms are ported here.

`sed_centered_diff_lhs` (oracle advance_microphys_module.F90:2188) is the implicit, centered-difference
mean-sedimentation operator for the d<hm>/dt equation:
    + (1/rho_ds) d( rho_ds <V_hm>|_t <hm>|_(t+1) ) / dz       (sign flipped on moving to the LHS).
<hm> lives on thermo (zt) levels, <V_hm> and rho_ds_zm on momentum (zm) levels; <hm> is interpolated to
zm with the zt2zm weights, multiplied by rho_ds_zm <V_hm>, and differenced over the central zt level.
The lower boundary sets <hm> at the surface (zm level 0) equal to the zt-level-0 value (no weight); the
upper boundary is no-flux (super-diagonal = 0).

The 3-band convention matches the solver / the budget code: band 0 = super-diagonal (couples hm[k+1]),
band 1 = main diagonal (hm[k]), band 2 = sub-diagonal (hm[k-1]) — i.e. Fortran kp1_tdiag/k_tdiag/
km1_tdiag = 1/2/3.

Validated bit-faithful-to-the-gate vs the rico `rrm_sd`/`Nrm_sd` budget stats (tests/test_kk_rico_oracle.py):
the stored budget term is `-sed_diff_lhs · <hm>` (3-point stencil with the advanced field).
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def sed_centered_diff_lhs(V_hm, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, w_above):
    """Mean-sedimentation implicit LHS (centered difference), all thermo levels at once.

    All inputs are (..., nz) arrays sharing a leading batch shape:
      V_hm, rho_ds_zm, w_above : MOMENTUM-level (nzm = nzt+1). V_hm is the sed velocity (<=0),
                                 rho_ds_zm the dry static density, w_above the zt2zm upper-thermo
                                 weight (= (zm[m]-zt[m-1])/(zt[m]-zt[m-1]); w_below = 1-w_above).
      invrs_rho_ds_zt, invrs_dzt : THERMO-level (nzt). invrs_dzt[k] = 1/(zm[k+1]-zm[k]).

    Returns (super, main, sub), each THERMO-level (..., nzt): the super-/main-/sub-diagonal bands.
    super[k] couples hm[k+1], main[k] hm[k], sub[k] hm[k-1]."""
    irzt, idzt = invrs_rho_ds_zt, invrs_dzt
    wa, wb = w_above, 1.0 - w_above
    nzt = irzt.shape[-1]
    c = irzt * idzt                                   # common prefactor, thermo (..., nzt)

    # Momentum arrays aligned to thermo level k: "k" -> index k, "kp1" -> index k+1.
    rzm_k, rzm_kp1 = rho_ds_zm[..., :nzt], rho_ds_zm[..., 1:nzt + 1]
    V_k, V_kp1 = V_hm[..., :nzt], V_hm[..., 1:nzt + 1]
    wa_k, wb_k = wa[..., :nzt], wb[..., :nzt]
    wa_kp1, wb_kp1 = wa[..., 1:nzt + 1], wb[..., 1:nzt + 1]

    # Interior (matches level 2..nzt-1 in Fortran 1-based; the bottom/top overwrite k=0, nzt-1).
    super = c * rzm_kp1 * V_kp1 * wa_kp1
    main = c * (rzm_kp1 * V_kp1 * wb_kp1 - rzm_k * V_k * wa_k)
    sub = -c * rzm_k * V_k * wb_k

    # Lower boundary (k=0): the -rho_ds_zm*V_hm term carries NO weight (hm_surface = hm[0]); sub=0.
    main = main.at[..., 0].set(c[..., 0] * (rzm_kp1[..., 0] * V_kp1[..., 0] * wb_kp1[..., 0]
                                            - rzm_k[..., 0] * V_k[..., 0]))
    sub = sub.at[..., 0].set(0.0)

    # Upper boundary (k=nzt-1): no flux through the model top -> super=0.
    super = super.at[..., nzt - 1].set(0.0)
    return super, main, sub


def term_turb_sed_lhs(Vhmphmp_impc, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, w_above):
    """Turbulent-sedimentation implicit LHS (centered difference). Oracle term_turb_sed_lhs
    (advance_microphys_module.F90:2683, l_upwind_diff_sed=.false. — the default/preferred branch).

    The semi-implicit turbulent-sed term is - (1/rho_ds) d( rho_ds <V_hm'h_m'> ) / dz with
    <V_hm'h_m'> = Vhmphmp_impc <h_m> + Vhmphmp_expc; the implicit part has the IDENTICAL flux-form
    discretization as the mean-sedimentation term (`sed_centered_diff_lhs`) with the momentum-level
    implicit covariance `Vhmphmp_impc` (= zt2zm of `KK_sed_vel_covars`' Vrrprrp_impc) playing the
    role of the velocity <V_hm>. Verified branch-by-branch identical to the Fortran. So this delegates
    to `sed_centered_diff_lhs`; the explicit part (Vhmphmp_expc) goes to the RHS.

    Returns (super, main, sub), the THERMO-level tridiagonal bands."""
    return sed_centered_diff_lhs(Vhmphmp_impc, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, w_above)


_C_K_HM = 0.75      # c_K_hm (parameters_tunable.F90:227)
_EPS = 1.0e-10      # constants_clubb eps = max(1e-10, machine eps)


def calculate_K_hm(wp2, Kh_zm, Skw_zm, hydromet, hydrometp2, hydromet_tol, gr,
                   c_K_hm=_C_K_HM):
    """Hydrometeor eddy-diffusivity K_hm on momentum levels (for the down-gradient turbulent-advection
    closure <w'hm'> = -K_hm d<hm>/dz). Oracle calculate_K_hm (advance_microphys_module.F90:3236),
    with l_use_non_local_diff_fac=.false. (the default — the K_gamma non-local factor is skipped):

      K_hm = c_K_hm Kh_zm (sqrt(hydrometp2) / max(zt2zm(hydromet), tol)) (1 + |Skw_zm|),

    then capped so |corr(w,hm)| <= 1: where |d<hm>/dz| > eps,
      K_hm = min( K_hm, sqrt(wp2) sqrt(hydrometp2) / |d<hm>/dz| ).
    K_hm = 0 at the lower/upper boundaries (not used). hydrometp2/Kh_zm/Skw_zm/wp2 are MOMENTUM-level;
    hydromet is THERMO-level (interpolated/differenced to momentum). Returns (..., nzm)."""
    from clubb_jax.src.CLUBB_core.grid_class import zt2zm, ddzt
    hm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet)  # momentum
    dhm = ddzt(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydromet)     # momentum: invrs_dzm[k]*(hm[k]-hm[k-1])
    sqrt_hmp2 = jnp.sqrt(jnp.maximum(hydrometp2, 0.0))
    K = (c_K_hm * Kh_zm * (sqrt_hmp2 / jnp.maximum(hm_zm, hydromet_tol))
         * (1.0 + jnp.abs(Skw_zm)))
    # Correlation cap (only where the gradient is non-negligible).
    adhm = jnp.abs(dhm)
    cap = jnp.sqrt(jnp.maximum(wp2, 0.0)) * sqrt_hmp2 / jnp.where(adhm > _EPS, adhm, 1.0)
    K = jnp.where(adhm > _EPS, jnp.minimum(K, cap), K)
    # Boundary conditions: K_hm = 0 at the bottom (k=0) and top (k=nzm-1).
    K = K.at[..., 0].set(0.0)
    K = K.at[..., -1].set(0.0)
    return K


def _turb_adv_lhs(K_hm, nu, rho_ds_zm, invrs_rho_ds_zt, gr):
    """Crank-Nicholson turbulent-advection LHS = (1/2) eddy diffusion, with the explicit k=0 lower-BC
    re-set (oracle :2096; byte-identical to diffusion_zt_lhs's own bottom row). Shared by
    microphys_lhs (the implicit half) and microphys_rhs (the explicit half). (3, ..., nzt)."""
    from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
    nu_arr = jnp.broadcast_to(jnp.asarray(nu, dtype=jnp.float64),
                              jnp.asarray(invrs_rho_ds_zt).shape[:-1])
    K_zt = jnp.zeros(nu_arr.shape + (gr.nzt,), dtype=jnp.float64)
    lhs_ta = 0.5 * diffusion_zt_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, gr,
        K_hm, K_zt, nu_arr, invrs_rho_ds_zt, rho_ds_zm,
    )
    idzt0 = gr.invrs_dzt[..., 0]; irzt0 = invrs_rho_ds_zt[..., 0]
    bc = 0.5 * irzt0 * (idzt0 * (K_hm[..., 1] + nu_arr) * rho_ds_zm[..., 1] * gr.invrs_dzm[..., 1])
    lhs_ta = lhs_ta.at[0, ..., 0].set(-bc)     # super
    lhs_ta = lhs_ta.at[1, ..., 0].set(bc)      # main
    lhs_ta = lhs_ta.at[2, ..., 0].set(0.0)     # sub
    return lhs_ta


def microphys_lhs(dt, K_hm, nu, wm_zt, V_hm, Vhmphmp_impc,
                  rho_ds_zm, invrs_rho_ds_zt, gr, w_above, l_sed=True):
    """Assemble the implicit LHS tridiagonal matrix for one hydrometeor's predictive equation.
    Oracle microphys_lhs (advance_microphys_module.F90:1564).

      lhs = 1/dt [main]
            + (1/2) diffusion_zt_lhs(K_hm, nu)       [Crank-Nicholson turbulent advection]
            + term_ma_zt_lhs(wm_zt)                  [mean advection]
            + sed_centered_diff_lhs(V_hm)        [mean sedimentation, if l_sed]
            + term_turb_sed_lhs(Vhmphmp_impc)    [turbulent sedimentation, if l_sed]

    Reuses the bit-faithful `diffusion_zt_lhs` (turb adv, the Fortran k=1 lower BC re-set is
    identical to its own bottom row, applied explicitly here for faithfulness) and
    `term_ma_zt_lhs` (upwind mean adv), plus the conservation-verified sedimentation operators.

    Inputs: K_hm, V_hm, Vhmphmp_impc, rho_ds_zm, w_above on MOMENTUM levels (..., nzm); wm_zt,
    invrs_rho_ds_zt on THERMO levels; nu the background diffusivity (..., or scalar broadcast);
    `gr` carries invrs_dzt/invrs_dzm. Returns (super, main, sub), each (..., nzt)."""
    from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs

    lhs_ta = _turb_adv_lhs(K_hm, nu, rho_ds_zm, invrs_rho_ds_zt, gr)
    lhs_ma = term_ma_zt_lhs(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        wm_zt,
        gr.weights_zt2zm,
        gr.invrs_dzt,
        gr.invrs_dzm,
        True,
        gr.grid_dir,
    )     # (3, ..., nzt) [super, main, sub]

    invrs_dzt = gr.invrs_dzt
    super = lhs_ta[0] + lhs_ma[0]
    main = lhs_ta[1] + lhs_ma[1] + 1.0 / dt
    sub = lhs_ta[2] + lhs_ma[2]

    if l_sed:
        ss, ms, subs = sed_centered_diff_lhs(V_hm, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, w_above)
        st, mt, subt = term_turb_sed_lhs(Vhmphmp_impc, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt, w_above)
        super = super + ss + st
        main = main + ms + mt
        sub = sub + subs + subt
    return super, main, sub


def lhs_budget_term(super, main, sub, hmm):
    """The budget tendency `-lhs · hmm` for a 3-band LHS (advance_microphys_module.F90:1518-1529):
      term[k] = -sub[k]*hmm[k-1] - main[k]*hmm[k] - super[k]*hmm[k+1]
    with hmm[k-1], hmm[k+1] clamped at the boundaries (km1=max(k-1,0), kp1=min(k+1,nzt-1))."""
    hmm = jnp.asarray(hmm)
    hmm_km1 = jnp.concatenate([hmm[..., :1], hmm[..., :-1]], axis=-1)      # hmm[max(k-1,0)]
    hmm_kp1 = jnp.concatenate([hmm[..., 1:], hmm[..., -1:]], axis=-1)      # hmm[min(k+1,nzt-1)]
    return -sub * hmm_km1 - main * hmm - super * hmm_kp1


def term_turb_sed_rhs(Vhmphmp_expc, rho_ds_zm, invrs_rho_ds_zt, invrs_dzt):
    """Explicit (RHS) part of the turbulent-sedimentation term (centered difference). Oracle
    term_turb_sed_rhs (advance_microphys_module.F90:3014, l_upwind_diff_sed=.false.):
      rhs[k] = - (1/rho_ds) d( rho_ds Vhmphmp_expc ) / dz
             = - invrs_rho_ds_zt[k] invrs_dzt[k] (F[k+1] - F[k]),   F[m] = rho_ds_zm[m] Vhmphmp_expc[m],
    a flux-divergence of the KNOWN explicit covariance (no <hm> factor); F is zeroed at the model top
    (momentum level nzm-1) for the no-flux upper boundary. Vhmphmp_expc is on MOMENTUM levels
    (= zt2zm of KK_sed_vel_covars' Vrrprrp_expc). Returns (..., nzt)."""
    nzt = jnp.asarray(invrs_rho_ds_zt).shape[-1]
    F = rho_ds_zm * Vhmphmp_expc
    F = F.at[..., nzt].set(0.0)                          # no flux through the model top
    return -invrs_rho_ds_zt * invrs_dzt * (F[..., 1:nzt + 1] - F[..., :nzt])


def microphys_rhs(dt, hmm, hmm_tndcy, K_hm, nu, Vhmphmp_expc,
                  rho_ds_zm, invrs_rho_ds_zt, gr, l_sed=True):
    """Assemble the explicit RHS vector for one hydrometeor's predictive equation. Oracle microphys_rhs
    (advance_microphys_module.F90:1912):
      rhs = hmm/dt + hmm_tndcy [microphysics source]
            - lhs_ta·hmm        [Crank-Nicholson EXPLICIT 1/2-diffusion of the current-step hmm]
            + term_turb_sed_rhs [explicit turbulent sedimentation, if l_sed]
    where lhs_ta = (1/2) diffusion_zt_lhs(K_hm,nu) (the SAME operator microphys_lhs uses implicitly).
    `hmm_tndcy` is the microphysics tendency (auto+accr+evap etc., from compute_kk_microphysics).
    Returns the RHS vector (..., nzt)."""
    lhs_ta = _turb_adv_lhs(K_hm, nu, rho_ds_zm, invrs_rho_ds_zt, gr)
    rhs = hmm / dt + hmm_tndcy + lhs_budget_term(lhs_ta[0], lhs_ta[1], lhs_ta[2], hmm)
    if l_sed:
        rhs = rhs + term_turb_sed_rhs(Vhmphmp_expc, rho_ds_zm, invrs_rho_ds_zt, gr.invrs_dzt)
    return rhs


def microphys_solve(lhs, rhs):
    """Solve the tridiagonal system for one hydrometeor mean field (advance_microphys_module.F90:microphys_solve).

    `lhs` is the (3, ..., nzt) [super, main, sub] band; returns the advanced <hm> (..., nzt). The Fortran's
    budget `stats_update`/`stats_finalize_budget` and error checks are not reproduced — the JAX path computes no
    microphys budget stats — mirroring the other extracted JAX solve routines (mfl_xm_solve, xp2_xpyp_solve)."""
    from clubb_jax.src.CLUBB_core.tridiag_lu_solver import tridiag_lu_solve
    return tridiag_lu_solve(rhs.shape[-1], lhs, rhs)


def advance_one_hydrometeor(dt, hmm, hmm_tndcy, K_hm, nu, wm_zt,
                            V_hm, Vhmphmp_impc, Vhmphmp_expc,
                            rho_ds_zm, invrs_rho_ds_zt, gr, w_above, l_sed=True):
    """Advance one hydrometeor mean field <hm> one timestep: assemble the implicit LHS + explicit RHS
    and solve the tridiagonal system. Oracle: microphys_lhs/microphys_rhs/microphys_solve
    (advance_microphys_module.F90). Returns the advanced <hm> (..., nzt).

    hmm is the current mean; hmm_tndcy the microphysics source (from compute_kk_microphysics). The sed
    velocities V_hm (mean) and Vhmphmp_impc/_expc (turbulent covariance impc/expc) are on momentum
    levels (zt2zm of KK_sedimentation / KK_sed_vel_covars). Hole-filling is applied separately."""
    sup, main, sub = microphys_lhs(dt, K_hm, nu, wm_zt, V_hm, Vhmphmp_impc,
                                   rho_ds_zm, invrs_rho_ds_zt, gr, w_above, l_sed)
    rhs = microphys_rhs(dt, hmm, hmm_tndcy, K_hm, nu, Vhmphmp_expc,
                        rho_ds_zm, invrs_rho_ds_zt, gr, l_sed)
    lhs = jnp.stack([sup, main, sub], axis=0)            # (3, ngrdcol, nzt) [super, main, sub]
    return microphys_solve(lhs, rhs)
