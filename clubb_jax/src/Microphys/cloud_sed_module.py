"""Cloud-droplet sedimentation — port of Microphys/cloud_sed_module.F90.

`cloud_drop_sed` accounts for gravitational settling of cloud droplets (Stokes
regime, lognormal size distribution; Ackerman et al. 2009, MWR, eq. 7). It is the
only "microphysics" some cases need (`l_cloud_sed=.true.` with
`microphys_scheme="none"`, e.g. atex_long, dycoms2_rf02_*). Produces the rcm/thlm
microphysics tendencies (`rcm_mc`, `thlm_mc`) added to the forcings next timestep.

Reference: clubb_release/src/Microphys/cloud_sed_module.F90:cloud_drop_sed
"""
from __future__ import annotations

import numpy as np
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.grid_class import zt2zm, ddzm
from clubb_jax.src.CLUBB_core.constants_clubb import rho_lw, Cp, Lv
from clubb_jax.src.CLUBB_core.tracer_numpy import _asarray  # REFACTOR B5: tracer-transparent np.asarray

_PI = float(np.pi)


def cloud_drop_sed(rcm, Ncm, rho_zm, rho, exner, sigma_g, gr):
    """Cloud-droplet sedimentation tendencies (cloud_drop_sed).

    Args:
        rcm:    (ngrdcol, nzt) mean cloud water mixing ratio [kg/kg].
        Ncm:    (ngrdcol, nzt) mean cloud droplet concentration [num/kg].
        rho_zm: (ngrdcol, nzm) density on momentum levels [kg/m^3].
        rho:    (ngrdcol, nzt) density on thermodynamic levels [kg/m^3].
        exner:  (ngrdcol, nzt) Exner function [-].
        sigma_g: geometric std dev of the droplet size distribution [-].
        gr:     grid object.

    Returns:
        (rcm_mc, thlm_mc, Fcsed): (ngrdcol, nzt) microphysics tendencies
        [kg/kg/s] and [K/s], plus the (ngrdcol, nzm) sedimentation flux [kg/(m^2 s)]
        for stats. (Fortran resets rcm_mc/thlm_mc to 0 before adding the
        sedimentation term, so for the cloud-sed-only path they equal the sed term.)
    """
    rcm_j    = jnp.asarray(rcm)
    Ncm_j    = jnp.asarray(Ncm)
    rho_zm_j = jnp.asarray(rho_zm)
    rho_j    = jnp.asarray(rho)
    exner_j  = jnp.asarray(exner)

    # Interpolate rcm, Ncm to momentum (zm) levels.
    rcm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, rcm_j)
    Ncm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm_j)

    # Sedimentation flux Fcsed on zm levels (positive downwards), only where there
    # is cloud (rcm>0 and Ncm>0). Use safe values inside the powers, then mask.
    valid    = (rcm_zm > 0.0) & (Ncm_zm > 0.0)
    Ncm_safe = jnp.where(Ncm_zm > 0.0, Ncm_zm, 1.0)
    rcm_safe = jnp.maximum(rcm_zm, 0.0)
    const    = 1.19e8 * jnp.exp(5.0 * (jnp.log(sigma_g)) ** 2)
    Fcsed_full = (const
                  * (3.0 / (4.0 * _PI * rho_lw * Ncm_safe * rho_zm_j)) ** (2.0 / 3.0)
                  * (rho_zm_j * rcm_safe) ** (5.0 / 3.0))
    Fcsed = jnp.where(valid, Fcsed_full, 0.0)
    # Boundary conditions: zero flux at the lowest and highest momentum levels.
    Fcsed = Fcsed.at[:, 0].set(0.0).at[:, -1].set(0.0)

    # d(rcm)/dt due to the sedimentation flux divergence (on zt levels).
    sed_rcm = (1.0 / rho_j) * ddzm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Fcsed)

    rcm_mc  = sed_rcm
    thlm_mc = -(Lv / (Cp * exner_j)) * sed_rcm
    # _asarray (REFACTOR B5): the body is already jnp (differentiable); only these returns severed the graph.
    return (_asarray(rcm_mc, dtype=np.float64),
            _asarray(thlm_mc, dtype=np.float64),
            _asarray(Fcsed, dtype=np.float64))
