"""Rayleigh-scattering optical depth + single-scattering albedo — port of rayle.F.

For shortwave band `ib`, the Rayleigh coefficient is a fixed `ri[ib]` except band 1 (ib=0 here),
where it is a quadratic in the cosine zenith angle u0. The optical depth per layer is
`tray = x · 29.267 · ppl · ln(pp[l+1]/pp[l])`; the single-scattering albedo `wray = 1`.
"""
import jax.numpy as jnp


def rayle(ib, amu0, ri, pp, ppl):
    """Rayleigh optical depth `tray` and SSA `wray` for shortwave band `ib` (0-based; Fortran ib is
    1-based, so band 1 ⇔ ib=0). `amu0` = (ncol,) cos(zenith); `ri` = (mbs,) band coefficients;
    `pp` = (ncol, nlm+1) level pressure [hPa]; `ppl` = (ncol, nlm) layer pressure [hPa].
    Returns (tray, wray), each (ncol, nlm)."""
    u0 = amu0
    if ib == 0:   # Fortran band 1
        x = -3.902860e-6 * u0 * u0 + 6.120070e-6 * u0 + 4.177440e-6
    else:
        x = jnp.full_like(u0, float(ri[ib]))
    x = x[:, None]                                       # (ncol, 1) broadcast over layers
    trp = 29.267 * ppl * jnp.log(pp[:, 1:] / pp[:, :-1])
    tray = x * trp
    wray = jnp.ones_like(tray)
    return tray, wray
