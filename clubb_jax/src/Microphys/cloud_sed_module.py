"""Cloud water sedimentation, mirroring Microphys/cloud_sed_module.F90.

JAX adaptation: column/level loops are array operations; inout values are returned; arrays carry
all columns.
"""

import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.grid_class import zt2zm, ddzm
from clubb_jax.src.CLUBB_core.constants_clubb import rho_lw, Cp, Lv, pi


# -----------------------------------------------------------------------------
def cloud_drop_sed(
    gr, ngrdcol, rcm, Ncm,        # In
    rho_zm, rho, exner, sigma_g,  # In
    stats, rcm_mc, thlm_mc,       # InOut
):
    """Account for cloud droplet sedimentation.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        rcm: Mean cloud water mixing ratio [kg/kg]
        Ncm: Mean cloud droplet concentration [num/kg]
        rho_zm: Density on momentum levels [kg/m^3]
        rho: Density on thermodynamic levels [kg/m^3]
        exner: Exner function [-]
        sigma_g: Geometric standard deviation of cloud droplets [-]
        stats: Immutable statistics state; return its updated value.
        rcm_mc: r_c tendency due to microphysics [kg/kg/s]
        thlm_mc: thlm tendency due to microphysics [K/s]
    """

    # Description:
    # Account for cloud droplet sedimentation.
    #
    # Sedimentation flux of cloud droplets should be treated by assuming a
    # log-normal size distribution of droplets falling in a Stoke's regime, in
    # which the sedimentation flux (Fcsed) is given by:
    #
    # Sedimentation Flux = constant
    #                      * [ ( 3 / ( 4 * pi * rho_lw * NcV ) )^(2/3) ]
    #                      * [ ( rho * rc )^(5/3) ]
    #                      * EXP[ 5 * ( ( LOG( sigma_g ) )^2 ) ];
    #
    # where constant = 1.19 x 10^8 (m^-1 s^-1) and sigma_g is the geometric
    # standard deviation of cloud droplets falling in a Stoke's regime.
    #
    # When written for a mass-dependent cloud droplet concentration, Nc:
    #
    # Sedimentation Flux = constant
    #                      * [ ( 3 / ( 4 * pi * rho_lw * Nc * rho ) )^(2/3) ]
    #                      * [ ( rho * rc )^(5/3) ]
    #                      * EXP[ 5 * ( ( LOG( sigma_g ) )^2 ) ].
    #
    # According to the above equation, sedimentation flux, Fcsed, is defined
    # positive downwards.  Therefore,
    #
    # (drc/dt)|_Fcsed = (1.0/rho) * d(Fcsed)/dz.
    # References:
    # http://journals.ametsoc.org/doi/abs/10.1175/2008MWR2582.1
    # -----------------------------------------------------------------------

    rcm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, rcm)
    Ncm_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ncm)
    # Define cloud water sedimentation flux on momentum levels.
    cloudy = (rcm_zm > 0.0) & (Ncm_zm > 0.0)
    # Mask inputs as well as outputs to keep inactive fractional powers finite.
    Fcsed = jnp.where(
        cloudy,
        1.19e8
        * (3.0 / (4.0 * pi * rho_lw * jnp.where(cloudy, Ncm_zm, 1.0) * rho_zm)) ** (2.0 / 3.0)
        * (rho_zm * jnp.maximum(rcm_zm, 0.0)) ** (5.0 / 3.0)
        * jnp.exp(5.0 * jnp.log(sigma_g) ** 2),
        0.0,
    )
    # Boundary conditions.
    Fcsed = Fcsed.at[:, 0].set(0.0).at[:, -1].set(0.0)
    # Find drc/dt due to cloud water sedimentation flux.
    sed_rcm = ddzm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Fcsed)
    sed_rcm = (1.0 / rho) * sed_rcm
    stats = stats.update("sed_rcm", sed_rcm)
    stats = stats.update("Fcsed", Fcsed)
    # + thlm/rtm_microphysics -- cloud water sedimentation.
    rcm_mc = rcm_mc + sed_rcm
    thlm_mc = thlm_mc - (Lv / (Cp * exner)) * sed_rcm
    return stats, rcm_mc, thlm_mc
