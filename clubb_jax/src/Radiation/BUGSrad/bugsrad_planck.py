"""Blackbody (Planck) band emission for BUGSrad longwave — port of bugsrad_planck.F90.

The first piece of the BUGSrad correlated-k radiation port (the prerequisite for gabls3 and the
~18 morrison+bugsrad cases). `planck` returns the blackbody flux emission [W/m^2] at the model
INTERFACES for one IR band, as a 5th-order polynomial fit in temperature (Horner form). There are
MBIR=12 IR bands, each with 6 fitted coefficients (`_B`).

Faithful to the Fortran (bugsrad_planck.F90:54-62):
  - interface 0 (top of model)  : poly(tt[0])
  - interfaces 1..NLM-1         : poly(0.5*(tt[k-1] + tt[k]))   (layer-mean temperature)
  - interface NLM (surface)     : poly(ts)
where tt is the NLM layer temperatures and bf has NLM+1 interfaces.
"""
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_MBIR = 12   # number of IR bands

# Fitted polynomial coefficients b(6, MBIR): _B[band] = the 6 coeffs (c0..c5) for that band.
# Listed band-by-band exactly as in the Fortran reshape((/.../), (/6, MBIR/)) (column-major →
# each source row is one band's 6 coefficients).
_B = np.array([
    [-25.889132,    0.75038381,   -0.87074567e-02,  0.50701144e-04, -0.14856755e-06,  0.17579587e-09],  # Band 1
    [ 25.397471,   -0.59596460,    0.53117737e-02, -0.21681758e-04,  0.36630792e-07, -0.11541419e-10],  # Band 2
    [ 57.891546,   -1.4745788,     0.14577775e-01, -0.68637478e-04,  0.14707480e-06, -0.98862337e-10],  # Band 3
    [ 21.837317,   -0.63194381,    0.71338812e-02, -0.38569394e-04,  0.95685257e-07, -0.76188561e-10],  # Band 4
    [  0.83155466, -0.15281669,    0.31020500e-02, -0.23768837e-04,  0.74605666e-07, -0.67494167e-10],  # Band 5
    [-19.432674,    0.37744942,   -0.22166529e-02,  0.11663914e-05,  0.22128830e-07, -0.28943829e-10],  # Band 6
    [-51.844021,    1.2280373,    -0.10600353e-01,  0.38135251e-04, -0.45111018e-07,  0.16679671e-10],  # Band 7
    [-31.210771,    0.85737498,   -0.87947387e-02,  0.39416747e-04, -0.67469797e-07,  0.43711306e-10],  # Band 8
    [ -5.4417604,   0.28970317,   -0.44571665e-02,  0.26395273e-04, -0.52111967e-07,  0.37627129e-10],  # Band 9
    [ 14.646543,   -0.25202253,    0.67234738e-03,  0.67552180e-05, -0.19815201e-07,  0.17221281e-10],  # Band 10
    [ 12.218584,   -0.31591213,    0.26032011e-02, -0.58878366e-05,  0.73276694e-08, -0.38798834e-11],  # Band 11
    [  1.0183416,  -0.79710154e-01, 0.13753393e-02, -0.40247214e-05,  0.63186167e-08, -0.41250652e-11],  # Band 12
], dtype=np.float64)


def _poly(T, c):
    """5th-order polynomial in T, Horner form: c0 + T(c1 + T(c2 + T(c3 + T(c4 + T·c5))))."""
    return c[0] + T * (c[1] + T * (c[2] + T * (c[3] + T * (c[4] + T * c[5]))))


def planck(tt, ts, nbir):
    """Blackbody emission bf [W/m^2] at the NLM+1 interfaces for IR band `nbir` (0-based, 0..11;
    Fortran nbir is 1-based). `tt` = (ncol, NLM) layer temperatures, `ts` = (ncol,) surface temp.
    Returns bf = (ncol, NLM+1)."""
    tt = jnp.asarray(tt, dtype=jnp.float64)
    ts = jnp.asarray(ts, dtype=jnp.float64)
    c = jnp.asarray(_B[nbir], dtype=jnp.float64)
    top = _poly(tt[:, :1], c)                              # interface 0: poly(tt[0])
    mid = _poly(0.5 * (tt[:, :-1] + tt[:, 1:]), c)         # interfaces 1..NLM-1: layer-mean temp
    sfc = _poly(ts[:, None], c)                            # interface NLM: poly(ts)
    return jnp.concatenate([top, mid, sfc], axis=1)
