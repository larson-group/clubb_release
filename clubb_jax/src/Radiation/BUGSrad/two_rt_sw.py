"""Two-stream shortwave radiative-transfer solver — port of two_rt_sw.F.

Computes the spectral direct/diffuse SW fluxes for band `ib` via delta-Eddington two-stream + the
adding technique, with the direct solar beam. The standard build sets `sel_rules_sw=.false.`
(bugs_rad.F:432), so only the FULL calculation is ported; `usenewexp` is not defined, so `exp()` is
the intrinsic exp (jnp.exp).

Per layer: delta-scale → direct-beam attenuation `exptau=exp(-τ*/μ0)` (cumulative product → `direct`),
the SW two-stream tr/rr, and the beam-driven layer sources sigu/sigd. Then the same top-down adding
(re/vd/vu/td) and bottom-up diffuse-flux propagation (fudif/fddif) as the LW solver, plus the direct
downward flux `fddir = μ0·slr·direct`.
"""
import jax.numpy as jnp
import jax.lax as lax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_EPS = 1.0e-02


def two_rt_sw(ib, slr, amu0, wc, asym, tau, asdir, asdif, mbs):
    """SW two-stream fluxes for 0-based band `ib` (0..mbs-1). `slr`/`amu0` = (ncol,) daylight
    fraction / cos(zenith); `wc`/`asym`/`tau` = (ncol, nlm); `asdir`/`asdif` = (ncol, mbs) direct/
    diffuse surface albedo. Returns (fudif, fddir, fddif), each (ncol, nlm+1)."""
    slr = jnp.asarray(slr, dtype=jnp.float64); amu0 = jnp.asarray(amu0, dtype=jnp.float64)
    wc = jnp.asarray(wc, dtype=jnp.float64); asym = jnp.asarray(asym, dtype=jnp.float64)
    tau = jnp.asarray(tau, dtype=jnp.float64)
    asdir = jnp.asarray(asdir, dtype=jnp.float64); asdif = jnp.asarray(asdif, dtype=jnp.float64)
    ncol, nlm = wc.shape
    mu0 = amu0[:, None]                                 # (ncol, 1)
    inv_mu0 = 1.0 / mu0

    # ---- per-layer SW delta-Eddington coefficients ----
    fact = asym * asym
    oms = ((1.0 - fact) * wc) / (1.0 - fact * wc)
    taus = (1.0 - fact * wc) * tau
    asy = asym / (1.0 + asym)
    exptau = jnp.exp(-taus * inv_mu0)                  # direct-beam transmission per layer
    t = 0.25 * (7.0 - oms * (4.0 + 3.0 * asy))
    r = -0.25 * (1.0 - oms * (4.0 - 3.0 * asy))
    kappa = jnp.sqrt(t ** 2 - r ** 2)
    rinf = r / (kappa + t)
    eggtau = jnp.exp(-kappa * taus)
    denom = 1.0 - rinf ** 2 * eggtau ** 2
    tr = (1.0 - rinf ** 2) * eggtau / denom
    rr = rinf * (1.0 - eggtau ** 2) / denom
    diff = kappa ** 2 - inv_mu0 ** 2
    near = jnp.abs(diff) < _EPS
    factb = jnp.where(near, 1.0 / _EPS, 1.0 / jnp.where(near, 1.0, diff))
    cc = oms * slr[:, None] * factb
    g3 = 0.5 - 0.75 * asy * mu0
    g4 = 1.0 - g3
    aa = g3 * (t - inv_mu0) + g4 * r
    bb = g4 * (t + inv_mu0) + g3 * r

    # direct beam: direct[0]=1, direct[l]=prod exptau[0..l-1]  (ncol, nlm+1)
    direct = jnp.concatenate([jnp.ones((ncol, 1)), jnp.cumprod(exptau, axis=1)], axis=1)
    dlay = direct[:, :-1]                               # direct at each layer's top (ncol, nlm)
    sigu = cc * ((aa - rr * bb) - aa * tr * exptau) * dlay
    sigd = cc * (-bb * tr + (bb - rr * aa) * exptau) * dlay

    # ---- 1. adding, top → down (carry re[l], vd[l]) ----
    def _add(carry, layer):
        re_l, vd_l = carry
        rr_l, tr_l, sigu_l, sigd_l = layer
        prop = 1.0 / (1.0 - re_l * rr_l)
        re_n = rr_l + tr_l ** 2 * re_l * prop
        vd_n = sigd_l + (tr_l * vd_l + tr_l * re_l * sigu_l) * prop
        vu_l = (rr_l * vd_l + sigu_l) * prop
        return (re_n, vd_n), (re_n, vd_n, vu_l, prop)

    z = jnp.zeros(ncol)
    (_, _), (re_nl, vd_nl, vu_l, td_l) = lax.scan(_add, (z, z), (rr.T, tr.T, sigu.T, sigd.T))
    re = jnp.concatenate([z[None], re_nl], axis=0).T   # (ncol, nlm+1)
    vd = jnp.concatenate([z[None], vd_nl], axis=0).T
    vu = vu_l.T; td = td_l.T

    # ---- 2. diffuse fluxes, bottom → up ----
    fudif_sfc = (asdif[:, ib] * vd[:, nlm] + asdir[:, ib] * slr * amu0 * direct[:, nlm]) \
        / (1.0 - asdif[:, ib] * re[:, nlm])

    def _flux(fu_l, idx):
        l = idx + 1
        fddif_l = re[:, l] * fu_l + vd[:, l]
        fu_lm1 = tr[:, idx] * fu_l * td[:, idx] + vu[:, idx]
        return fu_lm1, (fddif_l, fu_lm1)

    idxs = jnp.arange(nlm - 1, -1, -1)
    _, (fddif_desc, fudif_desc) = lax.scan(_flux, fudif_sfc, idxs)
    fddif = jnp.concatenate([jnp.zeros((ncol, 1)), fddif_desc[::-1].T], axis=1)        # fddif[0]=0
    fudif = jnp.concatenate([fudif_desc[::-1].T, fudif_sfc[:, None]], axis=1)

    # ---- 3. direct downward flux: fddir = μ0·slr·direct ----
    fddir = (amu0 * slr)[:, None] * direct
    return fudif, fddir, fddif
