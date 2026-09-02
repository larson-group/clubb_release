"""Cloud optical properties (ADT theory) for BUGSrad — port of cloudg.F.

`cloudg` returns the cloud optical depth `tcld`, single-scattering albedo `wcld`, and asymmetry
`asycld` for spectral band `ib`, from the cloud water/ice content `wcont` [g/m^3] and effective
radius `re` [µm], using Anomalous Diffraction Theory (Stephens et al. 1990 with Mitchell-1994
spherical corrections, modified-gamma m=0.5). Water clouds unless `flag` (ice).

REFACTOR A3 (iter8): the Fortran computes π in SINGLE precision (`acos(-1.)`) and truncates the
absorption core with `sngl(...)`. Those were deliberate single-precision artifacts the JAX once
reproduced for bit-faithfulness; under the numerical-accuracy standard we use float64 (simpler, ~1e-7
more accurate, well within Tier-C for gabls3 — the only bugsrad case).
  - band 1 (ib=0 here) computes extinction only (wcld=0.999999, asy=0.85, real refractive index).
"""
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_PI = float(np.pi)   # float64 π (REFACTOR A3: was the float32 value of π for bit-faithfulness)
_EPS = 1.0e-5


def cloudg(ib, pp, tt, wcont, re, pdist, cnrw, cniw, cnri, cnii, xlam, flag):
    """Cloud optical depth/SSA/asymmetry for 0-based band `ib`. `pp`=(ncol,nlm+1) level pres [hPa];
    `tt`/`wcont`/`re`=(ncol,nlm) temp [K]/content [g/m^3]/eff. radius [µm]; `pdist`=scalar dist. width;
    `cnrw/cniw/cnri/cnii/xlam`=(mb,) refractive indices + band centers; `flag`=ice if True.
    Returns (tcld, wcld, asycld), each (ncol, nlm)."""
    pp = jnp.asarray(pp, dtype=jnp.float64); tt = jnp.asarray(tt, dtype=jnp.float64)
    wcont = jnp.asarray(wcont, dtype=jnp.float64); re = jnp.asarray(re, dtype=jnp.float64)
    cnr = float(cnri[ib]) if flag else float(cnrw[ib])
    cni = float(cnii[ib]) if flag else float(cniw[ib])
    xl = float(xlam[ib])
    p0 = float(pdist); p1 = p0 + 1.0; p2 = p0 + 2.0
    f2 = p1 * p0; f3 = p2 * f2

    on = wcont > _EPS
    wcont_s = jnp.where(on, wcont, 1.0)                 # safe positive
    dz = 29.286 * jnp.log(pp[:, 1:] / pp[:, :-1]) * tt  # m
    rm = re / p2
    no = wcont_s / ((4.0 * _PI / 3.0) * f3 * 1.0e-6 * rm ** 3)
    area = 1.0e-6 * _PI * f2 * no * rm ** 2
    c0 = 2.0 * area; c1 = c0 / f2
    xm = 2.0 * _PI * rm / xl

    def _ext(um):
        # ext-core = real( p0/(um(um+1)^p1) + 1/(um^2 (um+1)^p0) - 1/um^2 ); um complex.
        return jnp.real(p0 / (um * (um + 1.0) ** p1)
                        + 1.0 / (um ** 2 * (um + 1.0) ** p0) - 1.0 / um ** 2)

    if ib == 0:   # Fortran band 1 — extinction only, real refractive index
        um = (2.0 * xm * (cnr - 1.0)).astype(jnp.complex128) * 1j
        ext = c0 + 2.0 * c1 * _ext(um)
        tcld = ext * dz
        wcld = jnp.full_like(tcld, 0.999999)
        asycld = jnp.full_like(tcld, 0.85)
    else:
        cm = complex(cnr, -cni)
        um = (2.0 * xm).astype(jnp.complex128) * ((cm - 1.0) * 1j)
        ext = c0 + 2.0 * c1 * _ext(um)
        vm = 4.0 * xm * cni                             # real
        expr = (p0 / (vm * (vm + 1.0) ** p1)
                + 1.0 / (vm ** 2 * (vm + 1.0) ** p0) - 1.0 / vm ** 2)
        abs_ = area + c1 * expr   # REFACTOR A3: float64 (was sngl(...) float32 truncation)
        tcld = ext * dz
        ext = jnp.where(ext < abs_, abs_, ext)
        wcld = (ext - abs_) / ext
        asycld = jnp.full_like(tcld, 0.85)

    tcld = jnp.where(on, tcld, 0.0)
    wcld = jnp.where(on, wcld, 0.0)
    asycld = jnp.where(on, asycld, 1.0)
    return tcld, wcld, asycld
