"""Correlated-k gaseous absorption helpers for BUGSrad — port of gases_ckd.F90 (algorithm only).

The correlated-k method (Fu 1991): per band × k-distribution interval `ig`, the absorption
coefficient is interpolated from precomputed coefficients `coefk(3, n_pressure)` in temperature
(`ln k = a + b·(T−T0) + c·(T−T0)²`) and linearly in pressure (`pscale` brackets the layer pressure
in the standard-pressure grid → ip1/ip2/pkd; `qk`/`qkio3` do the interpolation), then the gas
transmission optical depth is the absorber amount times k (`qoph2o`/`qopo3*`/`qopch4`/`qopn2o`/`qophc`).

This module ports the ALGORITHM helpers; the `gases` 18-band dispatch + the 1522-line `gases_ckd_data.h`
coefficient tables (hk*, c*h2o, c10ch4, …) are a separate (table-transcription) piece.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import (
    MOLAR_VOLUME, GRAVITY, MW_H2O, MW_O3, MW_DRY_AIR)

_CH4_CONC = 1.6e-6
_N2O_CONC = 0.28e-6


def pscale(ppl, stanp):
    """Bracket each layer pressure `ppl` (ncol, nlm) in the ascending standard-pressure grid `stanp`
    (n,). Returns (pkd, ip1, ip2): the linear weight and the 0-based bracketing indices (faithful to
    pscale's 3-case logic: below stanp[0] → ip1=ip2=0; ≥ stanp[-1] → top bracket; else interior)."""
    ppl = jnp.asarray(ppl, dtype=jnp.float64); stanp = jnp.asarray(stanp, dtype=jnp.float64)
    n = stanp.shape[0]
    j = jnp.clip(jnp.searchsorted(stanp, ppl, side='right') - 1, 0, n - 2)   # interior lower bracket
    below = ppl < stanp[0]
    above = ppl >= stanp[-1]
    pkd = jnp.where(below, ppl / stanp[0],
                    jnp.where(above, (ppl - stanp[n - 2]) / (stanp[n - 1] - stanp[n - 2]),
                              (ppl - stanp[j]) / (stanp[j + 1] - stanp[j])))
    ip1 = jnp.where(below, 0, jnp.where(above, n - 2, j))
    ip2 = jnp.where(below, 0, jnp.where(above, n - 1, j + 1))
    return pkd, ip1, ip2


def _qk_core(coefk, tt, pkd, ip1, ip2, t0, ylo, yhi, yhi_set):
    coefk = jnp.asarray(coefk, dtype=jnp.float64)         # (3, n)
    tt = jnp.asarray(tt, dtype=jnp.float64); pkd = jnp.asarray(pkd, dtype=jnp.float64)
    y1 = tt - t0
    y1 = jnp.where(y1 < ylo, ylo, y1)
    y1 = jnp.where(y1 > yhi, yhi_set, y1)                 # qkio3: >75 sets 70 (asymmetric, faithful)

    def kf(ip):
        a, b, c = coefk[0][ip], coefk[1][ip], coefk[2][ip]
        return jnp.exp(a + b * y1 + c * y1 ** 2)

    x1 = kf(ip1); x2 = kf(ip2)
    return jnp.where(ip1 == ip2, x1 * pkd, x1 + (x2 - x1) * pkd)


def qk(coefk, tt, pkd, ip1, ip2):
    """Interpolate the k-coefficient (Fu 1991): ln k = a + b·(T−245) + c·(T−245)² [clamp T−245 ∈
    [−65, 75]], then linear in pressure via (ip1, ip2, pkd)."""
    return _qk_core(coefk, tt, pkd, ip1, ip2, 245.0, -65.0, 75.0, 75.0)


def qkio3(coefk, tt, pkd, ip1, ip2):
    """k-coefficient for the IR O3 band: T−250 [clamp ∈ [−70, 75] but the >75 branch sets 70 — the
    Fortran's asymmetric clamp, replicated faithfully]."""
    return _qk_core(coefk, tt, pkd, ip1, ip2, 250.0, -70.0, 75.0, 70.0)


def qoph2o(fkg, dp, rmix):
    """H2O transmission optical depth (gases_ckd.F90:qoph2o)."""
    return jnp.asarray(fkg) * jnp.asarray(rmix) * jnp.asarray(dp) * MOLAR_VOLUME / MW_H2O * 10.0 / GRAVITY


def qopo3i(fkg, dp, o3mix):
    """O3 IR-band transmission (gases_ckd.F90:qopo3i)."""
    return jnp.asarray(fkg) * jnp.asarray(o3mix) * jnp.asarray(dp) * MOLAR_VOLUME / MW_O3 * 10.0 / GRAVITY


def qopo3s(fk, dp, o3mix):
    """O3 non-gray (SW) transmission with scalar k `fk` (gases_ckd.F90:qopo3s)."""
    return jnp.asarray(o3mix) * jnp.asarray(dp) * fk * MOLAR_VOLUME / MW_O3 * 10.0 / GRAVITY


def qophc(fkg, dp):
    """H2O/CO2 band-overlap transmission (gases_ckd.F90:qophc): tg = fkg·dp."""
    return jnp.asarray(fkg) * jnp.asarray(dp)


def qopch4(fkg, dp):
    """CH4 transmission (gases_ckd.F90:qopch4), ch4_conc = 1.6e-6 ppv (hardcoded in the Fortran)."""
    return jnp.asarray(fkg) * jnp.asarray(dp) * MOLAR_VOLUME * 10.0 / GRAVITY * _CH4_CONC / MW_DRY_AIR


def qopn2o(fkg, dp):
    """N2O transmission (gases_ckd.F90:qopn2o), n2o_conc = 0.28e-6 ppv (hardcoded in the Fortran)."""
    return jnp.asarray(fkg) * jnp.asarray(dp) * MOLAR_VOLUME * 10.0 / GRAVITY * _N2O_CONC / MW_DRY_AIR


# Per-band solar energy (Wm^-2) for the SW bands 1-6 that scale hk (gases_ckd.F90); LW bands use hk directly.
_SOLAR = {1: 619.618, 2: 484.295, 3: 149.845, 4: 48.7302, 5: 31.6576, 6: 5.79927}


def gases(ib, ig, dp, tt, rmix, o3mix, umco2, umch4, umn2o, pkd, ip1, ip2, pp):
    """Correlated-k gas transmission optical depth `tg` (ncol, nlm) + weight `hk` (scalar) for
    spectral band `ib` (1-based, 1..18) and k-interval `ig` (0-based). Dispatch over gases_ckd.F90's
    18 cases using the ported helpers + the parsed `gases_ckd_data` tables. `pkd/ip1/ip2` are the
    pressure-weighting from `pscale` (the band driver computes them with STANPS/STANPIR)."""
    from clubb_jax.src.Radiation.BUGSrad.gases_ckd_data import tables
    T = tables()

    def k3(name):     # qk on a 3D table's ig-slice
        return qk(T[name][:, :, ig], tt, pkd, ip1, ip2)

    def k2(name):     # qk on a 2D (ch4/n2o/hca/hcb) table
        return qk(T[name], tt, pkd, ip1, ip2)

    if ib == 1:
        tg = qopo3s(float(T['fk1o3'][ig]), dp, o3mix); hk = _SOLAR[1] * T['hk1'][ig]
    elif ib in (2, 3, 4, 5, 6):
        tg = qoph2o(k3(f'c{ib}h2o'), dp, rmix); hk = _SOLAR[ib] * T[f'hk{ib}'][ig]
    elif ib in (7, 8, 9, 13, 16, 17, 18):
        tg = qoph2o(k3(f'c{ib}h2o'), dp, rmix); hk = T[f'hk{ib}'][ig]
    elif ib in (10, 11):                                   # H2O + CH4 + N2O overlap
        tg1 = qoph2o(k3(f'c{ib}h2o'), dp, rmix)
        tg2 = qopch4(k2(f'c{ib}ch4'), dp)
        tg3 = qopn2o(k2(f'c{ib}n2o'), dp)
        tg = tg1 + tg2 / 1.6 * umch4 + tg3 / 0.28 * umn2o
        hk = T[f'hk{ib}'][ig]
    elif ib == 12:                                         # H2O + O3 overlap
        tg1 = qopo3i(qkio3(T['c12o3'][:, :, ig], tt, pkd, ip1, ip2), dp, o3mix)
        tg2 = qoph2o(k2('c12h2o'), dp, rmix)
        tg = tg1 + tg2; hk = T['hk12'][ig]
    elif ib in (14, 15):                                   # H2O + CO2 overlap (Fu approach 2)
        pq = jnp.where(jnp.asarray(pp)[:, :-1] >= 63.1, jnp.asarray(rmix), 0.0)   # pp level≥63.1 hPa
        fkga = k3(f'c{ib}hca'); fkgb = k3(f'c{ib}hcb')
        fkg = fkga / 330.0 * umco2 + pq * fkgb
        tg = qophc(fkg, dp); hk = T[f'hk{ib}'][ig]
    else:
        raise ValueError(f'gases: invalid band {ib}')
    return tg, float(hk)
