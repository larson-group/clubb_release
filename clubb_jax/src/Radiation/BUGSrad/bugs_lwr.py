"""Longwave band driver for BUGSrad — port of bugs_lwr.F (with -Dnooverlap, the CLUBB build).

Loops over the 12 LW spectral bands; per band: cloud optics (water + ice via `cloudg`, with the
asymmetry OVERRIDDEN by the fixed asym_wat/asym_ice tables where cloud is present), water-vapor
continuum (`gascon`), Planck emission (`planck`), gray-absorption combine (`comscp1`); then loops the
band's k-intervals: gas absorption (`gases`), non-gray combine (`comscp2`), clip τ≥0, two TWO-stream
solves (`two_rt_lw` for the cloudy and clear skies — the -Dnooverlap branch), accumulating
flux·hk into the all-sky/clear-sky up/down fluxes.
"""
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.Radiation.BUGSrad.cloudg import cloudg
from clubb_jax.src.Radiation.BUGSrad.gascon import gascon
from clubb_jax.src.Radiation.BUGSrad.bugsrad_planck import planck
from clubb_jax.src.Radiation.BUGSrad.comscp1 import comscp1
from clubb_jax.src.Radiation.BUGSrad.comscp2 import comscp2
from clubb_jax.src.Radiation.BUGSrad.gases_ckd import gases
from clubb_jax.src.Radiation.BUGSrad.two_rt_lw import two_rt_lw
from clubb_jax.src.Radiation.BUGSrad.gases_ckd_data import STANPIR, KG

_MB, _MBS, _MBIR = 18, 6, 12
_EPS, _TMAX, _PDIST = 1.0e-5, 340.0, 2.0
_KG_LW = KG[_MBS:]   # k-intervals for the 12 LW bands

# refractive-index + band-center tables (bugs_lwr.F:221-248); shared with bugs_swr.
CNRW = np.array([1.3422, 1.3281, 1.3174, 1.2901, 1.3348, 1.3700, 1.3191, 1.2821, 1.3160, 1.3030,
                 1.2739, 1.2319, 1.1526, 1.1981, 1.3542, 1.4917, 1.5463, 1.8718])
CNIW = np.array([6.4790e-9, 1.3417e-06, 1.2521e-4, 7.1533e-4, 4.2669e-2, 4.3785e-3, 1.3239e-2,
                 1.5536e-2, 5.3894e-2, 3.4346e-2, 3.7490e-2, 4.7442e-2, 1.2059e-1, 3.3546e-1,
                 4.1698e-1, 4.0674e-1, 3.6362e-1, 5.2930e-1])
CNRI = np.array([1.3266, 1.2986, 1.2826, 1.2556, 1.2963, 1.3956, 1.3324, 1.2960, 1.3121, 1.3126,
                 1.2903, 1.2295, 1.1803, 1.5224, 1.5572, 1.5198, 1.4993, 1.7026])
CNII = np.array([7.0696e-9, 9.1220e-7, 1.2189e-4, 5.7648e-4, 4.3144e-2, 8.2935e-3, 1.5540e-2,
                 2.5594e-2, 5.9424e-2, 5.1511e-2, 4.0325e-2, 4.7994e-2, 2.3834e-1, 3.0697e-1,
                 1.1852e-1, 4.3048e-2, 6.3218e-2, 1.5843e-1])
XLAM = np.array([0.45, 1.0, 1.6, 2.2, 3.0, 3.75, 4.878, 5.556, 6.452, 7.547, 8.511, 9.615,
                 11.236, 13.605, 16.529, 21.277, 29.412, 71.403])
ASYM_WAT = np.array([0.8200, 0.8547, 0.8619, 0.8683, 0.8723, 0.8703, 0.8566, 0.8040, 0.7463,
                     0.6579, 0.5103, 0.1279])
ASYM_ICE = np.array([0.8524, 0.8791, 0.9022, 0.8797, 0.8637, 0.8722, 0.8609, 0.8168, 0.7663,
                     0.6584, 0.6172, 0.3585])


def bugs_lwr(ts, ppl, dp, tt, rmix, cwrho, cirho, o3mix, pp, umco2, umch4, umn2o):
    """LW radiative fluxes. Inputs (ncol, nlm) unless noted: ts=(ncol,) surface temp; ppl/dp/tt/rmix/
    cwrho/cirho/o3mix layer fields; pp=(ncol, nlm+1) level pressure; umco2/umch4/umn2o gas concs.
    Returns (fdlw, fulw, fdlwcl, fulwcl), each (ncol, nlm+1) [down/up, all-sky/clear-sky]."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    ts, ppl, dp, tt, rmix = a(ts), a(ppl), a(dp), a(tt), a(rmix)
    cwrho, cirho, o3mix, pp = a(cwrho), a(cirho), a(o3mix), a(pp)
    ncol, nlm = tt.shape
    rew = jnp.full((ncol, nlm), 10.0); rei = jnp.full((ncol, nlm), 30.0)
    es = jnp.ones((ncol, _MBIR))                       # surface emissivity = 1
    from clubb_jax.src.Radiation.BUGSrad.gases_ckd import pscale
    pkd, ip1, ip2 = pscale(ppl, STANPIR)

    fdlw = jnp.zeros((ncol, nlm + 1)); fulw = jnp.zeros((ncol, nlm + 1))
    fdlwcl = jnp.zeros((ncol, nlm + 1)); fulwcl = jnp.zeros((ncol, nlm + 1))
    zero = jnp.zeros((ncol, nlm)); one = jnp.ones((ncol, nlm))

    for ib in range(_MBS + 1, _MB + 1):                # Fortran bands 7..18
        lw = ib - _MBS - 1                             # 0-based LW band index 0..11
        ib0 = ib - 1                                   # 0-based global band for cloudg/two_rt
        tcldw, wcldw, asycldw = cloudg(ib0, pp, tt, cwrho, rew, _PDIST, CNRW, CNIW, CNRI, CNII, XLAM, False)
        tcldi, wcldi, asycldi = cloudg(ib0, pp, tt, cirho, rei, _PDIST, CNRW, CNIW, CNRI, CNII, XLAM, True)
        asycldw = jnp.where(cwrho >= _EPS, ASYM_WAT[lw], asycldw)
        asycldi = jnp.where(cirho >= _EPS, ASYM_ICE[lw], asycldi)
        tgm = gascon(ib0, pp, ppl, dp, tt, rmix)
        bf = planck(tt, ts, lw)
        tau1, tauclr1, asym, asyclr, fwcld, fwclr = comscp1(
            zero, tcldi, tcldw, tgm, zero, zero, wcldi, wcldw, zero, one, asycldi, asycldw)
        for ig in range(_KG_LW[lw]):
            tg, hk = gases(ib, ig, dp, tt, rmix, o3mix, umco2, umch4, umn2o, pkd, ip1, ip2, pp)
            tau, tauclr, wc, wcclr = comscp2(tg, fwcld, fwclr, tau1, tauclr1)
            tau = jnp.maximum(tau, 0.0); tauclr = jnp.maximum(tauclr, 0.0)
            fug, fdg = two_rt_lw(ib0, wc, asym, tau, es, bf, _MBS)
            fugcl, fdgcl = two_rt_lw(ib0, wcclr, asyclr, tauclr, es, bf, _MBS)
            fdlw = fdlw + fdg * hk; fulw = fulw + fug * hk
            fdlwcl = fdlwcl + fdgcl * hk; fulwcl = fulwcl + fugcl * hk
    return fdlw, fulw, fdlwcl, fulwcl
