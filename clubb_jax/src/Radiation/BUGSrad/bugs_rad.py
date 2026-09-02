"""Main BUGSrad orchestration — port of bugs_rad.F (the -Dradoffline -Dnooverlap CLUBB build).

Does the cloud/gas preprocessing (bugs_rad.F:546-556), calls the LW band driver, rescales the solar
constant + calls the SW band driver for daytime columns (amu0≥0.01), and computes the LW/SW radiative
heating rates from the net-flux divergence. With -Dradoffline the radiation runs on the model grid
directly (nnp=nlm, no ghost top layer) and the TOA/surface diagnostic-flux block is skipped.

Preprocessing (radoffline ⇒ no ghost layer; nooverlap ⇒ cloud condensate weighted by cloud fraction):
  den   = ppl·100/(287·tl)        air density [kg/m³]   (ppl in mb)
  rmix  = ql/(1−ql)               vapour mixing ratio from specific humidity
  cwrho = den·1000·qcwl·acld      cloud-water content [g/m³], fraction-weighted  (likewise cirho)
The snow field qril and the overlap coefficients b1..b4 are UNUSED in the -Dradoffline -Dnooverlap
path (b1..b4 only feed the max/random-overlap two_rt_lw_iter variant), so they are omitted.

Heating rate (K/s): `rate[l] = -heat_fac·(Fnet[l] − Fnet[l+1])/dpl[l]`, `heat_fac = grav·0.01/cp`,
`Fnet = fup − fdn` (net upward flux), levels l/l+1 bracketing layer l, dpl in hPa (the 0.01 → Pa).
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.Radiation.BUGSrad.bugs_lwr import bugs_lwr
from clubb_jax.src.Radiation.BUGSrad.bugs_swr import bugs_swr
from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import GRAVITY as _GRAV, CP_DRY_AIR as _CP

_SOLAR_NORM = 1339.945   # the Σ of the BUGSrad SW band solar energies (bugs_rad.F:607 slr rescale)


def _heating(fup, fdn, fupcl, fdncl, dpl, heat_fac):
    fnet = fup - fdn                                   # net upward (ncol, nlm+1)
    fnetcl = fupcl - fdncl
    rate = -heat_fac * (fnet[:, :-1] - fnet[:, 1:]) / dpl       # (ncol, nlm)
    ratecl = -heat_fac * (fnetcl[:, :-1] - fnetcl[:, 1:]) / dpl
    return rate, ratecl


def bugs_rad(pl2, ppl, dpl, tl, ql, qcwl, qcil, o3l, ts, amu0, slr,
             alvdf, alndf, alvdr, alndr, acld, s0=1367.0, grav=_GRAV, cp=_CP,
             umco2=330.0, umch4=1.6, umn2o=0.28):
    """Compute the LW + SW radiative heating rates and fluxes (faithful bugs_rad.F argument order).
    `pl2`=(ncol,nlm+1) level pressure [mb]; `ppl`/`dpl`/`tl`/`ql`/`qcwl`/`qcil`/`o3l`/`acld`=(ncol,nlm)
    layer pressure [mb]/thickness [mb]/temp [K]/specific humidity/cloud-water mixing ratio/cloud-ice
    mixing ratio/ozone/cloud fraction; `ts`/`amu0`/`slr`/`alvdf`/`alndf`/`alvdr`/`alndr`=(ncol,) surface
    temp/cos zenith/daylight frac/diffuse-vis/diffuse-nIR/direct-vis/direct-nIR albedo; `s0`=solar
    constant. Returns dict with atl/asl/atlcl/aslcl (LW/SW all/clear heating rates, K/s) + the fluxes."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    pl2, ppl, dpl, tl = a(pl2), a(ppl), a(dpl), a(tl)
    ql, qcwl, qcil, o3l, acld = a(ql), a(qcwl), a(qcil), a(o3l), a(acld)
    ts, amu0, slr = a(ts), a(amu0), a(slr)
    ncol, nlm = tl.shape
    heat_fac = grav * 0.01 / cp

    # cloud/gas preprocessing (bugs_rad.F:546-556, radoffline ⇒ nnp=nlm, nooverlap ⇒ ×acld)
    den = ppl * 100.0 / (287.0 * tl)            # kg/m³ (ppl in mb → ×100 = Pa)
    rmix = ql / (1.0 - ql)                       # vapour mixing ratio
    cwrho = den * 1000.0 * qcwl * acld          # cloud-water content [g/m³], fraction-weighted
    cirho = den * 1000.0 * qcil * acld          # cloud-ice content   [g/m³]
    o3mix = o3l

    # LW (all columns)
    fdlw, fulw, fdlwcl, fulwcl = bugs_lwr(ts, ppl, dpl, tl, rmix, cwrho, cirho, o3mix, pl2,
                                          umco2, umch4, umn2o)
    atl, atlcl = _heating(fulw, fdlw, fulwcl, fdlwcl, dpl, heat_fac)

    # SW (daytime columns only: amu0 ≥ 0.01); slr rescaled to the actual solar constant
    slr_loc = slr * s0 / _SOLAR_NORM
    sw = bugs_swr(ts, amu0, slr_loc, alvdf, alndf, alvdr, alndr, ppl, dpl, tl, rmix,
                  cwrho, cirho, o3mix, pl2, umco2, umch4, umn2o)
    fdsw, fusw, fdswcl, fuswcl = sw[:4]
    day = (amu0 >= 0.01)[:, None]                       # zero SW at night
    fdsw = jnp.where(day, fdsw, 0.0); fusw = jnp.where(day, fusw, 0.0)
    fdswcl = jnp.where(day, fdswcl, 0.0); fuswcl = jnp.where(day, fuswcl, 0.0)
    asl, aslcl = _heating(fusw, fdsw, fuswcl, fdswcl, dpl, heat_fac)

    return dict(atl=atl, asl=asl, atlcl=atlcl, aslcl=aslcl,
                fdlw=fdlw, fulw=fulw, fdlwcl=fdlwcl, fulwcl=fulwcl,
                fdsw=fdsw, fusw=fusw, fdswcl=fdswcl, fuswcl=fuswcl,
                radvbc=sw[4], radvbccl=sw[5], radvdc=sw[6], radvdccl=sw[7],
                radnbc=sw[8], radnbccl=sw[9], radndc=sw[10], radndccl=sw[11])
