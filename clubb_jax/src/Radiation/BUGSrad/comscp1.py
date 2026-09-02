"""Combine the SCATTERING optical properties (Rayleigh + aerosol + ice/water cloud) into the
clear/all-sky band optical depth, the scattering-weighted SSA numerators `fwclr/fwcld`, and the
scattering-weighted asymmetry — port of comscp1.F.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def comscp1(taer, tcldi, tcldw, tgm, tray, waer, wcldi, wcldw, wray, asyaer, asycldi, asycldw):
    """Combine scattering optical properties (comscp1.F). All inputs (ncol, nlm). Returns
    (tccld1, tcclr1, asycld, asyclr, fwcld, fwclr)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    taer, tcldi, tcldw, tgm, tray = a(taer), a(tcldi), a(tcldw), a(tgm), a(tray)
    waer, wcldi, wcldw, wray = a(waer), a(wcldi), a(wcldw), a(wray)
    asyaer, asycldi, asycldw = a(asyaer), a(asycldi), a(asycldw)
    tcclr1 = tgm + tray + taer
    tccld1 = tcclr1 + tcldi + tcldw
    wwray = wray * tray; wwaer = waer * taer
    wwcldi = wcldi * tcldi; wwcldw = wcldw * tcldw
    fwclr = wwray + wwaer
    fwcld = fwclr + wwcldi + wwcldw
    asyclr = jnp.where(fwclr > 1.0e-10,
                       (asyaer * wwaer) / jnp.where(fwclr > 1.0e-10, fwclr, 1.0), 1.0)
    asycld = jnp.where(fwcld > 1.0e-10,
                       (asyaer * wwaer + asycldi * wwcldi + asycldw * wwcldw)
                       / jnp.where(fwcld > 1.0e-10, fwcld, 1.0), 1.0)
    return tccld1, tcclr1, asycld, asyclr, fwcld, fwclr
