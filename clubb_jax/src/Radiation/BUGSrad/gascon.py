"""Water-vapor continuum absorption (CKD2.4) for BUGSrad longwave — port of gascon.F.

`gascon` computes the H2O continuum optical depth `tgm` per layer for spectral band `ib`. SW bands
(global ib 1-6, here 0-5) have no continuum (`iflb=0`). For LW bands the continuum band is
`iflb[ib]` (12..1, reversed), and `parm_ckd24` is the CKD2.4 parameterized optical depth — a
regression `exp(c0 + c1·ln(amnt) + c2·T + c3·ln(patm) + c4·ph2o + c5·amnt + c6·ln(ph2o))` with a
small/large-amount region split (`ireg`) and a thin/thick-layer pathlength correction (`factor`).
`parm_ckd24` uses the INTRINSIC exp (NOT newexp — only the two-stream solvers use newexp).
"""
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.Radiation.BUGSrad.bugsrad_physconst import GRAVITY, R_D, R_STAR, MW_H2O, F_VIRT

# Fortran band (1-18) → CKD continuum band (1-12); 0 = no continuum (SW bands). 0-based array.
_IFLB = np.array([0, 0, 0, 0, 0, 0, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1])

# log-amount threshold for the small/large region split, per continuum band.
_H2OBND = np.array([-5, -3.5, -2.0, -2, -1., -4, -4, -4, -3, -3.5, -3, -2], dtype=np.float64)

# ck24_3(ncoef=7, nreg=2, nband=12) — here [band][reg][coef], transcribed band-by-band/region from gascon.F.
_CK24 = np.array([
    [[ 1.667e+00, 9.421e-01,-7.358e-03, 1.355e+00, 2.557e+03, 5.798e+01,-4.570e-01],   # band 1 reg 1
     [ 6.417e+00, 1.002e+00,-6.991e-03, 1.010e+00, 1.203e+01, 4.501e-02,-2.428e-02]],  # band 1 reg 2
    [[ 2.390e+00, 9.528e-01,-6.058e-03, 1.071e+00, 2.676e+02, 9.848e+00,-1.459e-01],   # band 2
     [ 4.849e+00, 1.002e+00,-6.910e-03, 8.961e-01, 1.635e+01, 2.115e-02, 7.243e-02]],
    [[ 2.326e+00, 9.720e-01,-6.551e-03, 8.739e-01, 6.984e+01, 8.346e-01, 4.824e-02],   # band 3
     [ 5.002e+00, 1.005e+00,-9.286e-03, 6.222e-01, 1.168e+01, 3.611e-03, 3.148e-01]],
    [[-4.865e+00, 8.455e-01,-6.911e-03, 1.475e+00, 2.905e+02, 7.078e+00,-6.846e-01],   # band 4
     [ 4.596e+00, 1.012e+00,-1.152e-02, 5.713e-01, 1.270e+01,-1.395e-03, 3.447e-01]],
    [[-5.396e+00, 8.596e-01,-8.479e-03, 1.619e+00, 1.664e+02, 3.236e+00,-7.782e-01],   # band 5
     [ 7.478e+00, 1.007e+00,-1.963e-02, 2.771e-01, 6.021e+00,-4.489e-03, 6.709e-01]],
    [[ 1.262e+00, 2.347e-01,-2.360e-02, 1.655e-01, 5.068e+02, 2.462e+01, 3.920e-01],   # band 6
     [ 9.334e+00, 1.002e+00,-2.429e-02, 3.575e-02, 2.751e-01,-1.189e-03, 9.593e-01]],
    [[-1.222e+00, 5.423e-01,-2.327e-02, 5.197e-01, 6.423e+02, 5.038e+01, 1.502e-01],   # band 7
     [ 8.506e+00, 1.000e+00,-2.339e-02, 8.891e-03,-6.805e-01,-1.639e-04, 9.917e-01]],
    [[-3.638e+00, 8.534e-01,-1.344e-02, 6.816e-01, 5.385e+02, 4.428e+01,-6.366e-03],   # band 8
     [ 6.921e+00, 1.002e+00,-1.974e-02, 6.350e-02, 6.838e-01,-1.121e-03, 9.237e-01]],
    [[-2.329e+00, 7.893e-01,-2.588e-03, 1.017e+00, 1.525e+02, 1.029e+01,-1.486e-01],   # band 9
     [ 6.742e-01, 1.008e+00,-3.376e-03, 9.105e-01, 1.074e+01,-3.307e-03, 5.741e-02]],
    [[-1.677e+00, 9.173e-01,-5.780e-03, 1.504e+00, 7.886e+02, 2.288e+01,-5.999e-01],   # band 10
     [ 3.396e+00, 1.005e+00,-3.433e-03, 1.012e+00, 7.635e+00, 3.010e-03,-2.418e-02]],
    [[ 7.943e-01, 9.260e-01,-5.050e-03, 1.141e+00, 2.221e+02, 1.021e+01,-2.246e-01],   # band 11
     [ 3.356e+00, 1.002e+00,-4.719e-03, 9.578e-01, 6.164e+00, 1.186e-03, 2.264e-02]],
    [[-5.874e+00, 7.060e-01,-1.532e-03, 1.141e+00, 1.463e+02, 6.534e+00,-4.308e-01],   # band 12
     [ 4.709e-01, 1.010e+00,-6.067e-03, 8.513e-01, 1.161e+01,-6.629e-03, 8.885e-02]],
], dtype=np.float64)


def parm_ckd24(iband0, amnt, patm, temp, dz):
    """CKD2.4 continuum optical depth for continuum band `iband0` (0-based, 0..11). amnt [g/cm^2],
    patm [atm], temp [K], dz [km] — same shape; returns the optical depth. Vectorized."""
    amnt = jnp.asarray(amnt, dtype=jnp.float64); patm = jnp.asarray(patm, dtype=jnp.float64)
    temp = jnp.asarray(temp, dtype=jnp.float64); dz = jnp.asarray(dz, dtype=jnp.float64)
    c = jnp.asarray(_CK24[iband0], dtype=jnp.float64)   # (2, 7) [reg][coef]
    h2b = float(_H2OBND[iband0])
    # thin/thick-layer pathlength correction
    factor = jnp.where(dz < 0.25, 0.25 / dz, jnp.where(dz > 1.50, 1.50 / dz, 1.0))
    dz1 = jnp.where(dz < 0.25, 0.25, jnp.where(dz > 1.50, 1.50, dz))
    amnt1 = amnt * factor
    ph2o = amnt1 * (R_STAR * 1.e4 * temp) / (dz1 * 1.0e5 * MW_H2O * 1.01325e6)
    ireg = (jnp.log(amnt1) > h2b)                       # False→region 1 (idx 0), True→region 2 (idx 1)
    creg = jnp.where(ireg[..., None], c[1], c[0])       # (..., 7)
    patmx = jnp.log(patm)
    tau_log = (creg[..., 0] + creg[..., 1] * jnp.log(amnt1) + creg[..., 2] * temp
               + creg[..., 3] * patmx + creg[..., 4] * ph2o + creg[..., 5] * amnt1
               + creg[..., 6] * jnp.log(ph2o))
    return jnp.exp(tau_log) / factor


def gascon(ib, pp, ppl, dp, tt, rmix):
    """H2O continuum optical depth `tgm` (ncol, nlm) for 0-based global band `ib` (0..17). 0 where
    the band has no continuum (SW) or rmix<=0. `pp`=(ncol,nlm+1) level pres [hPa], `ppl`/`dp`/`tt`/
    `rmix`=(ncol,nlm) layer pres/thickness [hPa]/temp [K]/vapor mixing ratio [kg/kg]."""
    pp = jnp.asarray(pp, dtype=jnp.float64); ppl = jnp.asarray(ppl, dtype=jnp.float64)
    dp = jnp.asarray(dp, dtype=jnp.float64); tt = jnp.asarray(tt, dtype=jnp.float64)
    rmix = jnp.asarray(rmix, dtype=jnp.float64)
    cb = int(_IFLB[ib])
    if cb == 0:
        return jnp.zeros_like(ppl)
    iband0 = cb - 1
    on = rmix > 0.0
    rmix_s = jnp.where(on, rmix, 1.0e-10)               # safe positive so the logs never see <=0
    amnt = 10.0 * dp * rmix_s / GRAVITY                 # g/cm^2
    patm = ppl / 1013.25                                # atm
    tv = tt * (1.0 + F_VIRT * rmix_s)
    dz = (R_D / GRAVITY) * tv * jnp.log(pp[:, 1:] / pp[:, :-1]) * 0.001   # km
    tgm = parm_ckd24(iband0, amnt, patm, tt, dz)
    return jnp.where(on, tgm, 0.0)
