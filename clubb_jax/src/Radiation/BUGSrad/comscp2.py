"""Add the gas absorption optical depth `tg` to the scattering optical depth (from comscp1) and
form the final clear/all-sky single-scattering albedos (capped at 0.999999) — port of comscp2.F.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def comscp2(tg, fwcld, fwclr, tccld1, tcclr1):
    """Add gas absorption `tg` and form the final SSAs (comscp2.F). Returns (tccld, tcclr, wccld,
    wcclr); SSAs capped at 0.999999."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    tg, fwcld, fwclr, tccld1, tcclr1 = a(tg), a(fwcld), a(fwclr), a(tccld1), a(tcclr1)
    tcclr = tcclr1 + tg
    tccld = tccld1 + tg
    wcclr = jnp.where(tcclr > 0.0, fwclr / jnp.where(tcclr > 0.0, tcclr, 1.0), 0.0)
    wcclr = jnp.minimum(0.999999, wcclr)
    wccld = jnp.where(tccld > 0.0, fwcld / jnp.where(tccld > 0.0, tccld, 1.0), 0.0)
    wccld = jnp.minimum(0.999999, wccld)
    return tccld, tcclr, wccld, wcclr
