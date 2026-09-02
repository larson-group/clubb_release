"""Two-stream longwave radiative-transfer solver — port of two_rt_lw.F.

Computes the spectral up/down LW fluxes for band `ib` using the delta-Eddington two-stream method
with the adding technique. The standard BUGSrad build sets `sel_rules_lw=.false.` (bugs_rad.F:432),
so the selection-rules emission-only branch is dead — only the FULL calculation is ported. The build
does NOT define `usenewexp` (CPPDEFS has no -Dusenewexp), so `exp()` is the INTRINSIC exp (jnp.exp),
NOT the newexp approximation.

Per layer: delta-scale (fact=asym²) → diffusivity-two-stream tr/rr + the layer sources sigu/sigd.
Then two sequential recursions (lax.scan): the top-down adding (re/vd/vu/td) and the bottom-up flux
propagation (fu/fd).
"""
import jax.numpy as jnp
import jax.lax as lax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_DIFFAC = 2.0


def two_rt_lw(ib, wc, asym, tau, es, bf, mbs):
    """LW two-stream fluxes for 0-based band `ib`. `wc`/`asym`/`tau` = (ncol, nlm) SSA/asymmetry/
    optical depth; `es` = (ncol, mbir) surface emissivity; `bf` = (ncol, nlm+1) Planck at levels;
    `mbs` = number of SW bands (ibms = ib - mbs indexes es). Returns (fu, fd), each (ncol, nlm+1)."""
    wc = jnp.asarray(wc, dtype=jnp.float64); asym = jnp.asarray(asym, dtype=jnp.float64)
    tau = jnp.asarray(tau, dtype=jnp.float64); bf = jnp.asarray(bf, dtype=jnp.float64)
    es = jnp.asarray(es, dtype=jnp.float64)
    ncol, nlm = wc.shape
    ibms = ib - mbs
    bfl, bfr = bf[:, :-1], bf[:, 1:]                    # bf[l], bf[l+1]   (ncol, nlm)

    # ---- per-layer delta-Eddington two-stream coefficients ----
    fact = asym * asym
    oms = ((1.0 - fact) * wc) / (1.0 - fact * wc)
    taus = (1.0 - fact * wc) * tau
    beta0 = (4.0 + asym) / (8.0 * (1.0 + asym))
    t = _DIFFAC * (1.0 - oms * (1.0 - beta0))
    r = _DIFFAC * oms * beta0
    kappa = jnp.sqrt(t ** 2 - r ** 2)
    rinf = r / (kappa + t)
    ggtau = kappa * taus
    eggtau = jnp.exp(-ggtau)
    denom = 1.0 - rinf ** 2 * eggtau ** 2
    tr = (1.0 - rinf ** 2) * eggtau / denom
    rr = rinf * (1.0 - eggtau ** 2) / denom
    # layer sources sigu/sigd: thin-layer limit vs full
    thin = taus < 0.1481e-2
    sigu_thin = 0.5 * _DIFFAC * (bfl + bfr) * taus
    aa = (t + r) * (1.0 - rr) - (1.0 + rr - tr) / jnp.where(thin, 1.0, taus)
    bb = -(t + r) * tr + (1.0 + rr - tr) / jnp.where(thin, 1.0, taus)
    cc = _DIFFAC * (1.0 - oms) / kappa ** 2
    sigu = jnp.where(thin, sigu_thin, cc * (aa * bfl + bb * bfr))
    sigd = jnp.where(thin, sigu_thin, cc * (bb * bfl + aa * bfr))

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
    layers = (rr.T, tr.T, sigu.T, sigd.T)              # (nlm, ncol) each → scan over layer axis
    (_, _), (re_nl, vd_nl, vu_l, td_l) = lax.scan(_add, (z, z), layers)
    # re/vd have nlm+1 entries (re[0]=vd[0]=0, then re[1..nlm] from the scan)
    re = jnp.concatenate([z[None], re_nl], axis=0).T   # (ncol, nlm+1)
    vd = jnp.concatenate([z[None], vd_nl], axis=0).T
    vu = vu_l.T                                        # (ncol, nlm)
    td = td_l.T

    # ---- 2. fluxes, bottom → up (carry fu[l]) ----
    fu_sfc = es[:, ibms] * bf[:, nlm]                  # fu at the surface level (index nlm)
    # iterate l = nlm .. 1: fd[l]=re[l]fu[l]+vd[l]; fu[l-1]=tr[l-1] fu[l] td[l-1] + vu[l-1]
    def _flux(fu_l, idx):
        # idx = l-1 (0-based layer index for tr/td/vu), l = idx+1
        l = idx + 1
        fd_l = re[:, l] * fu_l + vd[:, l]
        fu_lm1 = tr[:, idx] * fu_l * td[:, idx] + vu[:, idx]
        return fu_lm1, (fd_l, fu_lm1)

    idxs = jnp.arange(nlm - 1, -1, -1)                 # nlm-1 .. 0  (l = nlm .. 1)
    _, (fd_desc, fu_desc) = lax.scan(_flux, fu_sfc, idxs)
    # fd_desc[k] = fd at level nlm-k → reversed gives levels 1..nlm; fu_desc[k] = fu at level nlm-1-k.
    fd_top = jnp.zeros((ncol, 1))                       # fd[0] = 0 (top boundary)
    fd = jnp.concatenate([fd_top, fd_desc[::-1].T], axis=1)           # (ncol, nlm+1)
    fu = jnp.concatenate([fu_desc[::-1].T, fu_sfc[:, None]], axis=1)  # (ncol, nlm+1)
    return fu, fd
