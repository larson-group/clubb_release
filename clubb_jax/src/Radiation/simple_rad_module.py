r"""JAX port of ``Radiation/simple_rad_module.F90``.

The simplified LW radiation code used by both various idealized cases is now
consolidated into the simple_rad subroutine below.  The comments here are
in regard to the differences between the COAMPS-LES simulations done in the
past and how the idealized radiation is handled here in CLUBB.
-dschanen 29 Sep 2010

SPECIAL METHOD USED TO CALCULATE RADIATION
Grid descriptions by Adam Smith, 27 June 2006

In order to verify our CLUBB simulations are working properly, we
have first developed a series of 3D simulations using the COAMPS-LES
model.  This large-eddy simulation (LES) simulation uses specific
methods to calculate radiation, subsidence, and other microphysical
processes.  To make the two models simluate clouds as closely as
possible, we use the same radiation scheme in both models.

In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement
all radiation code.  This allows the subroutine to be duplicated
exactly in many different models.  However, the subroutine uses the
COAMPS vertical grid.  Therefore, for CLUBB to implement this code
correctly, we must modify some of our variable profiles before
calling the radiation subroutine.

The following diagram describes the differences in model grids
(see "ADDITIONAL NOTE" for important grid information):

      COAMPS-LES                                   CLUBB

 k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)

 k= 1   (t) ----------                   ---------- k=kk        (t)

 k= 2   (m) ----------                   ---------- k=kk        (m)

 k= 2   (t) ----------                   ---------- k=kk-1      (t)

           .                  .                     .
           .                  .                     .
           .                  .                     .

 k=kk-1 (m) ----------  m = momentum     ---------- k=3         (m)
                                 levels
 k=kk-1 (t) ----------  t = thermo       ---------- k=2         (t)
                                 levels
 k=kk   (m) ----------                   ---------- k=2         (m)

 k=kk   (t) ----------  kk = number of   ---------- k=1         (t)
                             vertical
 k=kk+1 (m) ----------       heights     ---------- k=1         (m)

//////////////////////// MODEL SURFACE /////////////////////////////

ADDITIONAL NOTE: In order to reconcile the COAMPS and CLUBB grids,
                 the uppermost level of the COAMPS grid is ignored,
                 eliminating the need to add level kk+1 to the CLUBB
                 grid.  For the purposes of this code, the COAMPS
                 grid levels are re-indexed to start from level 1
                 at the uppermost useful COAMPS grid level (which
                 was previously referred to as level 2 in the above
                 diagram and is at the same altitude as CLUBB
                 level kk).  Likewise, the lowermost COAMPS grid
                 level is now indexed at kk (rather than kk+1 in the
                 above diagram).  Brian Griffin; May 10, 2008.

Also, the COAMPS grid indices are numbered from the top of the model
downward, while the CLUBB grid indices are numbered from the bottom
up.  Therefore, since we are using a COAMPS radiation scheme, we
flip moisture and temperature profiles that are passed into the
rad_lwsw subroutine.  The rad scheme will produce results in using
the COAMPS grid scheme, so all radiation output will be flipped
back to the CLUBB grid before being applied to the model.

Finally, since the COAMPS scheme does not have a gridpoint below
model surface, we add that point to all radiative output files once
they are converted back to CLUBB setup.  This allows all averages and
calculations to be done correctly.


Computation of radiative fluxes on staggered grid
Comments by Michael Falk, 16 February 2005.

Frad (and its components Frad_LW and Frad_SW) should be computed on
w points, not on mass points, which is apparent from its formulation
and from its location in stats_sw instead of stats_sm.  The grid
looks like this:


-----Frad----------------------------------    k = 1  (w level)
    /    \            |-dwm
-LWP------radht----------------------------    k = 1  (mass level)
    \    /            |
-----Frad----------------------------------    k = 2  (w level)
    /    \
-LWP------radht----------------------------    k = 2  (mass level)
    \    /
-----Frad----------------------------------    k = 3  (w level)
    /    \
-LWP------radht----------------------------    k = 3  (mass level)

If you consider Frad to take place on mass levels, then computing
LWP is a forward difference and is only first-order accurate, while
if Frad computed in between LWP levels, it is a centered difference
which is second-order accurate.

The coding implementation requires that Frad depend on LWP(k) and
LWP(k-1) since the w level for a given k is at a higher altitude
than the mass level.  radht, back on mass levels, depends on Frad(k)
and Frad(k+1).

ADDITIONAL NOTE: For clarification of terminology, a w level on the
                 COAMPS grid is equivalent to a momentum (m) level
                 on the CLUBB grid, and a mass level on the COAMPS
                 grid is equivalent to a thermodynamic (t) level on
                 the CLUBB grid.  Brian Griffin; May 10, 2008.

Additionally, these computations assume that the distance between
mass levels (dsigma) is constant, and that the w levels (spaced by
dsigmw) always fall exactly halfway in between the mass levels.  If
this is not the case, consider dwm to be the distance between a w
level and the mass level below it, and dmw to be the distance
between a mass level and the w level below it.  Then, the
formulation for Frad_LW, for instance, would use a weighted average:

(dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
which, for dwm always == dmw, reduces to
(1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
which is identical to the current formulation.
((lwp(k)+lwp(k-1))/2)

ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
                 compatible with the use of a stretched
                 (or unevenly-spaced) grid, as well as with the use
                 of an evenly-spaced grid.  Interpolation functions
                 are used to compute any weighted averages, rather
                 than using general numbers such as (1/2), which is
                 compatible only with an evenly-spaced grid.
                 Brian Griffin; May 10, 2008.

JAX adaptation:
    ``RadiationParameters`` is the immutable equivalent of the source
    ``parameters_radiation`` module variables. Fortran output and inout values
    are returned explicitly, in source ownership order. LBA file I/O remains a
    host initialization boundary; all per-timestep calculations are JAX.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, eps
from clubb_jax.src.CLUBB_core.interpolation import (
    lin_interpolate_two_points,
    zlinterp_fnc,
)

configure_jax_precision()

# These variables are for the LBA radiation
# Constant Parameters
lba_ntimes = 36
lba_nzrad = 33
# Altitudes        [m]
# Radiative tendencies     [K/s]


@partial(jax.jit, static_argnames=("ngrdcol",))
def simple_rad(
    gr, ngrdcol,
    rho,  # Density on thermodynamic grid  [kg/m^3]
    rho_zm,  # Density on momentum grid       [kg/m^3]
    rtm,  # Total water mixing ratio       [kg/kg]
    rcm,  # Cloud water mixing ratio       [kg/kg]
    exner,  # Exner function.                [-]
    # JAX adaptation: immutable source ``parameters_radiation`` module state.
    stats, radiation_parameters,
):
    """Description:
        A simplified radiation driver

    References:
        None

    Returns ``stats`` (inout), ``Frad_LW`` (out), and ``radht_LW`` (out), in
    the order of the Fortran argument list after its input fields.
    """
    # External

    # Constant parameters
    ls_div = 3.75e-6

    # Input Variables

    # Output Variables

    # Local Variables

    # ---- Begin Code ----

    LWP = liq_water_path(ngrdcol, gr.nzm, gr.nzt, rho, rcm, gr.invrs_dzt)

    if radiation_parameters.F1 > eps:
        Frad_LW = (
            radiation_parameters.F0 * jnp.exp(-radiation_parameters.kappa * LWP)
            + radiation_parameters.F1
            * jnp.exp(-radiation_parameters.kappa * (LWP[:, 0:1] - LWP))
        )
    else:  # Mathematically equivalent to the above, but computationally cheaper
        Frad_LW = radiation_parameters.F0 * jnp.exp(-radiation_parameters.kappa * LWP)

    if radiation_parameters.l_rad_above_cloud:
        # Find the height of the isotherm rtm = 8.0 g/kg.
        above_threshold = rtm > 8.0e-3
        k_iso = jnp.sum(
            jnp.cumprod(above_threshold.astype(jnp.int32), axis=1), axis=1
        )
        # JAX gathers require an in-bounds index. The source assumes this
        # isotherm lies strictly inside the column; keep that precondition and
        # only clip the index used to express its valid-data calculation.
        k_iso = jnp.clip(k_iso, 1, gr.nzt - 1)

        k_iso_minus_one = k_iso - 1
        rtm_high = jnp.take_along_axis(rtm, k_iso[:, None], axis=1)[:, 0]
        rtm_low = jnp.take_along_axis(rtm, k_iso_minus_one[:, None], axis=1)[:, 0]
        zt_high = jnp.take_along_axis(gr.zt, k_iso[:, None], axis=1)[:, 0]
        zt_low = jnp.take_along_axis(gr.zt, k_iso_minus_one[:, None], axis=1)[:, 0]
        z_i = lin_interpolate_two_points(
            8.0e-3, rtm_high, rtm_low, zt_high, zt_low,
        )

        # Compute the Heaviside step function for z - z_i.
        dz = gr.zm - z_i[:, None]
        Heaviside = jnp.where(dz < -eps, 0.0, jnp.where(dz > eps, 1.0, 0.5))

        # ``jnp.where`` evaluates both branches, unlike the source ``if``.
        # Clamp only the fractional-power input to preserve the source's
        # positive-Heaviside calculation without creating negative fractional
        # powers for levels below the inversion.
        dz_positive = jnp.maximum(dz, 0.0)
        Frad_LW = Frad_LW + (
            rho_zm * Cp * ls_div * Heaviside
            * (0.25 * dz_positive ** (4.0 / 3.0)
               + z_i[:, None] * dz_positive ** (1.0 / 3.0))
        )

        if stats.l_sample:
            # Update inversion-height statistics used by surface-radiation diagnostics.
            stats = stats.update("z_inversion", z_i)

    # Compute the radiative heating rate.
    # The radiative heating rate is defined on thermodynamic levels.
    radht_LW = (
        (1.0 / exner) * (-1.0 / (Cp * rho))
        * (Frad_LW[:, 1:] - Frad_LW[:, :-1]) * gr.invrs_dzt
    )

    return stats, Frad_LW, radht_LW


@partial(jax.jit, static_argnames=("ngrdcol",))
def simple_rad_bomex(gr, ngrdcol):
    """Description:
        Compute radiation as per the GCSS BOMEX specification.

    References:
        <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>

    Returns the Fortran ``radht`` output array.
    """
    # Output Variables

    # Local Variables

    # ---- Begin Code ----

    # Radiative theta-l tendencysimple_rad_bomex
    return jnp.where(
        (gr.zt >= 0.0) & (gr.zt < 1500.0),
        -2.315e-5,
        jnp.where(
            (gr.zt >= 1500.0) & (gr.zt < 2500.0),
            # From bomex specification, section 3.4
            -2.315e-5 + 2.315e-5 * (gr.zt - 1500.0) / (2500.0 - 1500.0),  # Known magic number
            0.0,
        ),
    )


@partial(jax.jit, static_argnames=("ngrdcol",))
def simple_rad_lba(
    gr, ngrdcol,
    time_current,  # Current time of model run   [s]
    time_initial,  # Start time of model run     [s]
    # JAX adaptation: source ``lba_zrad``/``lba_krad`` module state.
    radiation_parameters,
):
    """Description:
        Compute radiation For the LBA TRMM case.  Uses a prescribed formula and
        interpolates with respect to time.

    References:
        None

    ``radiation_parameters`` supplies the source module's ``lba_zrad`` and
    ``lba_krad`` arrays; the returned value is Fortran ``radht``.
    """
    # Input Variables

    # Output Variables

    # Local Variables

    time = time_current - time_initial

    # Calculate radiative heating rate
    radhtz = lba_radhtz(time, radiation_parameters.lba_krad)

    # Radiative theta-l tendency
    return jax.vmap(
        lambda zt: zlinterp_fnc(zt, radiation_parameters.lba_zrad, radhtz)
    )(gr.zt)


# JAX adaptation: source ``iunit`` is a Fortran file-unit argument; direct host
# path I/O has no equivalent unit and retains only the source data location.
def simple_rad_lba_init(file_path: str | Path):
    """Initialize the module by reading the LBA forcing data.

    Description:
        This subroutine initializes the module by reading in forcing
        data used in the tndcy subroutine.

    This host I/O boundary replaces the source ``iunit`` plus
    ``file_read_1d``/``file_read_2d`` interface and returns ``lba_zrad`` then
    ``lba_krad``, the source module variables it initializes.
    """
    # Constant Parameters

    # ---- Begin Code ----
    path = Path(file_path)
    lba_zrad = np.fromstring((path / "lba_heights.dat").read_text().replace(",", " "), sep=" ")
    lba_krad = np.fromstring((path / "lba_rad.dat").read_text().replace(",", " "), sep=" ")
    if lba_zrad.size != lba_nzrad:
        raise ValueError(f"lba_heights.dat: expected {lba_nzrad}, got {lba_zrad.size}")
    if lba_krad.size != lba_nzrad * lba_ntimes:
        raise ValueError(
            f"lba_rad.dat: expected {lba_nzrad * lba_ntimes}, got {lba_krad.size}"
        )
    return jnp.asarray(lba_zrad), jnp.asarray(lba_krad.reshape(lba_nzrad, lba_ntimes))


@partial(jax.jit, static_argnames=("ngrdcol", "nzm", "nzt"))
def liq_water_path(
    ngrdcol,  # Number of grid columns in the model
    nzm,  # Number of momentum vertical levels in the model
    nzt,  # Number of thermodynamic vertical levels in the model
    rho,  # Air Density                      [kg/m^3]
    rcm,  # Cloud water mixing ratio         [kg/kg]
    invrs_dzt,  # Inverse of distance per level    [1/m]
):
    """Description:
        Compute liquid water path

    References:
        None

    Returns the source function result ``liq_water_path`` on momentum levels.
    """
    # Input Variables

    # Output Variables

    # ---- Begin Code ----

    # Liquid water path is defined on the intermediate model levels between the
    # rcm and rho levels (i.e. the momentum levels in CLUBB).
    contribution = rcm * rho / invrs_dzt

    def accumulate(LWP_above, contribution_at_level):
        LWP = LWP_above + contribution_at_level
        return LWP, LWP

    _, LWP_reversed = jax.lax.scan(
        accumulate,
        jnp.zeros((ngrdcol,), dtype=rho.dtype),
        contribution[:, ::-1].T,
    )
    LWP = LWP_reversed.T[:, ::-1]
    return jnp.concatenate((LWP, jnp.zeros((ngrdcol, 1), dtype=rho.dtype)), axis=1)


# TODO(port-mirror): ``simple_rad_lba`` has this time-interpolation loop inline
# in Fortran. JAX needs a separately jitted callable to represent its dynamic
# branch without host scalar conversion; keep it local to this module.
@jax.jit
def lba_radhtz(time, lba_krad):
    """Calculate the LBA table value at elapsed model time."""
    i1 = jnp.clip(jnp.floor(time / 600.0).astype(jnp.int32), 1, lba_ntimes - 1)
    i2 = i1 + 1
    a = (time - 600.0 * i1) / (600.0 * i2 - 600.0 * i1)
    interpolated = a * (lba_krad[:, i2 - 1] - lba_krad[:, i1 - 1]) + lba_krad[:, i1 - 1]
    return jnp.where(
        time <= 600.0,
        lba_krad[:, 0],
        jnp.where(time >= lba_ntimes * 600.0, lba_krad[:, lba_ntimes - 1], interpolated),
    )


__all__ = [
    "simple_rad", "simple_rad_bomex", "simple_rad_lba", "simple_rad_lba_init",
    "liq_water_path",
]
