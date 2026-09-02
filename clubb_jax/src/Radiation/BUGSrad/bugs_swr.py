"""Shortwave band driver for BUGSrad — port of bugs_swr.F (with -Dnooverlap, the CLUBB build).

Parallels bugs_lwr but for the 6 SW bands: Rayleigh scattering (`rayle`), cloud optics (`cloudg`,
asymmetry overridden by the SW asym_wat/asym_ice tables), gray combine (`comscp1`); then per
k-interval: gas absorption (`gases`, using the temperature-capped `ttem=min(tmax,tt)`), non-gray
combine (`comscp2`), two `two_rt_sw` solves (cloudy + clear, with the direct beam + surface albedos),
accumulating direct+diffuse fluxes and the surface visible/near-IR net radiation.
"""
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.Radiation.BUGSrad.rayle import rayle
from clubb_jax.src.Radiation.BUGSrad.cloudg import cloudg
from clubb_jax.src.Radiation.BUGSrad.comscp1 import comscp1
from clubb_jax.src.Radiation.BUGSrad.comscp2 import comscp2
from clubb_jax.src.Radiation.BUGSrad.gases_ckd import gases, pscale
from clubb_jax.src.Radiation.BUGSrad.two_rt_sw import two_rt_sw
from clubb_jax.src.Radiation.BUGSrad.gases_ckd_data import STANPS, KG
from clubb_jax.src.Radiation.BUGSrad.bugs_lwr import CNRW, CNIW, CNRI, CNII, XLAM   # shared tables

_MBS = 6
_EPS, _TMAX, _PDIST = 1.0e-5, 340.0, 2.0
_KG_SW = KG[:_MBS]
_RI = np.array([0.9022e-5, 0.5282e-6, 0.5722e-7, 0.1433e-7, 0.4526e-8, 0.1529e-8])
ASYM_WAT_SW = np.array([0.8625, 0.8469, 0.8287, 0.8182, 0.9472, 0.7630])
ASYM_ICE_SW = np.array([0.8678, 0.8640, 0.8653, 0.8615, 0.9526, 0.8293])


def bugs_swr(ts, amu0, slr, alvdf, alndf, alvdr, alndr, ppl, dp, tt, rmix,
             cwrho, cirho, o3mix, pp, umco2, umch4, umn2o):
    """SW radiative fluxes + surface net radiation. (ncol,) scalars: ts, amu0 (cos zenith), slr
    (daylight frac), alvdf/alndf/alvdr/alndr (visible/near-IR diffuse/direct surface albedo). (ncol,
    nlm): ppl/dp/tt/rmix/cwrho/cirho/o3mix; pp=(ncol,nlm+1). Returns (fdsw, fusw, fdswcl, fuswcl,
    radvbc, radvbccl, radvdc, radvdccl, radnbc, radnbccl, radndc, radndccl)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    ts, amu0, slr = a(ts), a(amu0), a(slr)
    ppl, dp, tt, rmix = a(ppl), a(dp), a(tt), a(rmix)
    cwrho, cirho, o3mix, pp = a(cwrho), a(cirho), a(o3mix), a(pp)
    ncol, nlm = tt.shape
    rew = jnp.full((ncol, nlm), 10.0); rei = jnp.full((ncol, nlm), 30.0)
    ttem = jnp.minimum(_TMAX, tt)
    # surface albedos: band 1 = visible, bands 2-6 = near-IR
    asdir = jnp.stack([a(alvdr)] + [a(alndr)] * 5, axis=1)   # (ncol, 6)
    asdif = jnp.stack([a(alvdf)] + [a(alndf)] * 5, axis=1)
    pkd, ip1, ip2 = pscale(ppl, STANPS)

    fdsw = jnp.zeros((ncol, nlm + 1)); fusw = jnp.zeros((ncol, nlm + 1))
    fdswcl = jnp.zeros((ncol, nlm + 1)); fuswcl = jnp.zeros((ncol, nlm + 1))
    z = jnp.zeros((ncol, nlm)); one = jnp.ones((ncol, nlm))
    rad = {k: jnp.zeros(ncol) for k in ('vbc', 'vbccl', 'vdc', 'vdccl', 'nbc', 'nbccl', 'ndc', 'ndccl')}

    for ib in range(1, _MBS + 1):                       # Fortran SW bands 1..6
        ib0 = ib - 1
        tray, wray = rayle(ib0, amu0, _RI, pp, ppl)
        tcldw, wcldw, asycldw = cloudg(ib0, pp, tt, cwrho, rew, _PDIST, CNRW, CNIW, CNRI, CNII, XLAM, False)
        tcldi, wcldi, asycldi = cloudg(ib0, pp, tt, cirho, rei, _PDIST, CNRW, CNIW, CNRI, CNII, XLAM, True)
        asycldw = jnp.where(cwrho >= _EPS, ASYM_WAT_SW[ib0], asycldw)
        asycldi = jnp.where(cirho >= _EPS, ASYM_ICE_SW[ib0], asycldi)
        tau1, tauclr1, asym, asyclr, fwcld, fwclr = comscp1(
            z, tcldi, tcldw, z, tray, z, wcldi, wcldw, wray, one, asycldi, asycldw)
        for ig in range(_KG_SW[ib0]):
            tg, hk = gases(ib, ig, dp, ttem, rmix, o3mix, umco2, umch4, umn2o, pkd, ip1, ip2, pp)
            tau, tauclr, wc, wcclr = comscp2(tg, fwcld, fwclr, tau1, tauclr1)
            fugdif, fdgdir, fdgdif = two_rt_sw(ib0, slr, amu0, wc, asym, tau, asdir, asdif, _MBS)
            fugcldif, fdgcldir, fdgcldif = two_rt_sw(ib0, slr, amu0, wcclr, asyclr, tauclr, asdir, asdif, _MBS)
            fdsw = fdsw + (fdgdir + fdgdif) * hk; fusw = fusw + fugdif * hk
            fdswcl = fdswcl + (fdgcldir + fdgcldif) * hk; fuswcl = fuswcl + fugcldif * hk
            if ib == 1:                                 # visible
                rad['vbc'] = rad['vbc'] + fdgdir[:, nlm] * hk
                rad['vbccl'] = rad['vbccl'] + fdgcldir[:, nlm] * hk
                rad['vdc'] = rad['vdc'] + fdgdif[:, nlm] * hk
                rad['vdccl'] = rad['vdccl'] + fdgcldif[:, nlm] * hk
            else:                                       # near-IR (bands 2-6)
                rad['nbc'] = rad['nbc'] + fdgdir[:, nlm] * hk
                rad['nbccl'] = rad['nbccl'] + fdgcldir[:, nlm] * hk
                rad['ndc'] = rad['ndc'] + fdgdif[:, nlm] * hk
                rad['ndccl'] = rad['ndccl'] + fdgcldif[:, nlm] * hk
    return (fdsw, fusw, fdswcl, fuswcl, rad['vbc'], rad['vbccl'], rad['vdc'], rad['vdccl'],
            rad['nbc'], rad['nbccl'], rad['ndc'], rad['ndccl'])
