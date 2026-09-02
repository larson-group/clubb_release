"""JAX port of ``Radiation/rad_lwsw_module.F90`` shortwave radiation.

REFERENCES (for rad_lwsw and sunray_sw):
Bluestein, H. B., 1992: Synoptic-Dynamic Meteorology in
      Midlatitudes, Volume I.  Oxford Univerity Press, 431 pp.
Duynkerke, P.G. et al., 2005: Observations and numerical simulations
      of the diurnal cycle of the EUROCS stratocumulus case.  In
      press, special EUROCS issue, Quart. Jour. Roy. Met. Soc.
Salby, M. L., 1996: Fundamentals of Atmospheric Physics.  Academic
      Press, 627 pp.
Shettle, E. P., and J. A. Weinman, 1970: The transfer of solar
      irradiance through inhomogenous turbid atmospheres evaluated
      by Eddington's approximation.  J. Atm. Sci., 27, 1048-1055.
Stevens, B. et al., 2005: Evaluation of large-eddy simulations via
      observations of nocturnal marine stratocumulus.  Submitted to
      Mon. Wea. Rev.
Wallace, J. M., and P. V. Hobbs, 1977: Atmospheric Science: An
      Introductory Survey.  Academic Press, 467 pp.
"""

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import rho_lw, three_halves
from clubb_jax.src.CLUBB_core.interpolation import lin_interpolate_two_points

configure_jax_precision()


@partial(jax.jit, static_argnames=("ngrdcol", "nzt", "l_center"))
def sunray_sw(
    ngrdcol, nzt,
    rcm, rho, xi_abs, dzt,
    zm,  # Altitude of momentum levels (CLUBB bottom-up)      [m]
    zt,  # Altitude of thermodynamic levels (CLUBB bottom-up) [m]
    radius, A, gc, Fs0, omega, l_center,
):
    """Compute shortwave radiative flux.

    Description:
    for CLEX altocumulus case
    Written by Geert Lenderink for implementation of Shettle and
    Weinman's formulation for radiative flux into the EUROCS
    stratocumulus case.

    Adapted by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
    others for COAMPS and CLUBB.

    Subroutine to compute shortwave radiative flux.

    The code for sunray_sw was slightly reconstructed in order to more
    closely follow Geert Lenderink's code for the Duynkerke et al.
    EUROCS case.  The original formulation used in that paper comes from
    Shettle and Weinman case.

    Tau is now computed inside this routine, as it's not needed in
    nov11_rad.  Tau is also computed for each layer instead of as the
    total optical depth from the top of the domain to the bottom of the
    layer being computed.  This makes tau zero everywhere outside of
    cloud.  F_diff and F_dir are also zero outside of cloud.

    Comments by Michael Falk, 26 January 2005.

    ADDITION TO COMMENT BY ADAM SMITH, 28 July 2006.
    We are attempting to compare CLUBB results with the COAMPS-LES model.
    To keep the two codes as similar as possible, this subroutine is a
    near-duplicate of the rad_lwsw subroutine  used in COAMPS_LES.  I
    have  modified variable declarations as needed to match CLUBB
    standards, and sections that do not apply to CLUBB have been removed.

    For a full explanation of how this subroutine is implemented in CLUBB,
    please see the nov11_altocu_tndcy or jun25_altocu_tndcy subroutines
    within gcss.F.

    This version operates natively on CLUBB's bottom-up grid ordering
    (k=1 at surface, k=nzt at top for thermo levels, k=nzm=nzt+1 at
    top for momentum levels) and handles multiple grid columns.
    In CLUBB's grid, thermo level k sits between momentum levels k and
    k+1.  Optical depth (taupath) accumulates from the top of the
    domain downward, i.e. from momentum level nzm down to level 1.

    REFERENCES:
    see subroutine rad_lwsw.

    Returns the Fortran ``Frad_SW`` output on momentum levels, with shape
    ``(ngrdcol, nzt + 1)``.
    """
    # Input variables

    # Output variables

    # Local Variables

    # -----------------------------------------------------------------------
    # CONSTANTS/PARAMETERS
    #
    # values added by Michael Falk and Adam Smith
    #
    # ff     : gc^2, denoted "g^2" in Duynkerke.            Unit: none
    # gcde   : Delta-Eddington transformation of gc.
    #          Notated g-prime in Duynkerke eqn.18.         Unit: None
    #
    # -----------------------------------------------------------------------
    ff = gc * gc
    gcde = gc / (1.0 + gc)

    # -----------------------------------------------------------------------
    #
    # Computation of tau and omega variables
    #
    #
    # tauc   : column total optical depth                   Unit: none
    # taupath: column total Delta-Eddington optical depth.  Unit: none
    # omega  : single-scattering albedo                     Unit: none
    # omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none
    # taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none
    # taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none
    #
    # Comments by Michael Falk, 26 January 2005
    #
    # -----------------------------------------------------------------------
    omegade = (1.0 - ff) * omega / (1.0 - omega * ff)

    # Compute per-layer optical depth and column total.
    # tauc accumulates over k per column, so parallelize over i only.
    tau = three_halves * rcm * rho * dzt / radius / rho_lw  # Optical depth of an incremental layer.        [-]
    tauc = jnp.sum(tau, axis=1)  # Column total optical depth                    [-]
    taude = (1.0 - omega * ff) * tau  # Delta-Eddington transformation of tau.        [-]

    # -----------------------------------------------------------------------
    #
    # Computation of constants for radiative transfer equations
    #
    # These variables come from Duynkerke eqn.20, which, with slight
    # modifications, matches Shettle and Weinman between eqns.12b and 13.
    # Duynkerke uses slightly different formulations than Shettle and
    # Weinman:
    #
    # F0(S&W)    = F0(Duynkerke)*pi.
    # alpha(S&W) = alpha(Duynkerke)*F0(S&W).
    # beta(S&W)  = beta(Duynkerke)*F0(S&W).
    # c1(S&W)    = c1(Duynkerke)*F0(S&W).
    # c2(S&W)    = c2(Duynkerke)*F0(S&W).
    # x3(S&W)    = x3(Duynkerke)*F0(S&W).
    # y3(S&W)    = y3(Duynkerke)*F0(S&W).
    #
    # F0 is divided out of each term in several equations, and then
    # reintroduced when the actual radiative flux is computed.  The
    # computations here follow Lenderink's original sunray_sw code which
    # uses the Duynkerke formulations for these variables.
    #
    # x1     : term 1 in k equation                         Unit: none
    # x2     : term 2 in k equation                         Unit: none
    # rk     : k equation                                   Unit: none
    # rk2    : k equation squared                           Unit: none
    # x3     : term in denominator of alpha and beta        Unit: none
    # rp     : p equation                                   Unit: none
    # alpha  : alpha equation                               Unit: none
    # beta   : beta equation                                Unit: none
    #
    # The following variables are used by Lenderink to solve for
    # Duynkerke's parameters C1 and C2.  They are all dimensionless.
    #
    # rtt    : 2/3.
    # exmu0  : exponential term used in S&W eqn.14-- originally from
    #          eqn 1 (also Salby 9.35) in the source function for SW
    #          radiation
    # expk   : one of the coefficients of C2 on the left side of Shettle
    #          and Weinman eqn.14, originally from eqn.12a and eqn.12b
    # exmk   : reciprocal of expk, one of the coefficients of C1 on the
    #          left side of Shettle and Weinman eqn.14, originally from
    #          eqn.12a and eqn.12b
    # xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13
    # xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13
    # ap23b  : right side of Shettle and Weinman eqn.13
    # t1     : the other coefficient of C1 on left side of Shettle and
    #          Weinman eqn.14
    # t2     : the other coefficient of C2 on left side of Shettle and
    #          Weinman eqn.14
    # t3     : the coefficient of exmu0 on the right side of Shettle and
    #          Weinman eqn.14
    #
    # Comments by Michael Falk, 26 January 2005
    #
    # -----------------------------------------------------------------------
    x1 = 1.0 - omegade * gcde
    x2 = 1.0 - omegade
    rk = jnp.sqrt(3.0 * x2 * x1)
    xi_abs2 = xi_abs * xi_abs
    rk2 = rk * rk
    x3 = 4.0 * (1.0 - rk2 * xi_abs2)
    rp = jnp.sqrt(3.0 * x2 / x1)
    alpha = 3.0 * omegade * xi_abs2 * (1.0 + gcde * x2) / x3
    beta = 3.0 * omegade * xi_abs * (1.0 + 3.0 * gcde * xi_abs2 * x2) / x3

    rtt = 2.0 / 3.0
    xp23p = 1.0 + rtt * rp
    xm23p = 1.0 - rtt * rp
    ap23b = alpha + rtt * beta
    t1 = 1.0 - A - rtt * (1.0 + A) * rp
    t2 = 1.0 - A + rtt * (1.0 + A) * rp
    t3 = (1.0 - A) * alpha - rtt * (1.0 + A) * beta + A * xi_abs

    # -----------------------------------------------------------------------
    #
    # Shettle and Weinman 13 and 14, adapted into Duynkerke, give two
    # equations and two unknowns which can be solved by linear combination
    # (C2) and then substitution (C1).
    #
    # exmu0, expk, exmk, taucde are per-column since they depend on the
    # column total optical depth (tauc).
    #
    # Comments by Michael Falk, 26 January 2005
    #
    # -----------------------------------------------------------------------
    taucde = (1.0 - omega * ff) * tauc  # D-E column total optical depth                [-]
    exmu0 = jnp.exp(-taucde / xi_abs)  # exp(-taucde/xi_abs)                           [-]
    expk = jnp.exp(rk * taucde)  # exp(rk*taucde)                                [-]
    exmk = 1.0 / expk  # 1/expk                                        [-]
    c2 = (xp23p * t3 * exmu0 - t1 * ap23b * exmk) / (
        xp23p * t2 * expk - xm23p * t1 * exmk
    )
    c1 = (ap23b - c2 * xm23p) / xp23p  # Shettle-Weinman coefficients                  [-]

    # -----------------------------------------------------------------------
    #
    # Computation of diffuse and direct shortwave flux
    #
    # F_diff is the first term in Duynkerke eqn.19.  The F0 and pi
    # constants which are divided out in the Duynkerke formulation are
    # reintroduced here.
    #
    # Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8,
    # where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,
    # where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in
    # Shettle and Weinman eqn.12b should actually be a lowercase p.
    #
    # F_dir is the second term in Duynkerke eqn.19.
    #
    # The negative sign for F_diff and F_dir is related to the definition
    # of which direction is a positive flux.
    #
    # Comments by Michael Falk, 26 January 2005
    #
    # -----------------------------------------------------------------------

    # -----------------------------------------------------------------------
    #
    # Computation of shortwave fluxes on staggered grid
    #
    # For a full explanation see the "Computation of radiative fluxes on
    # staggered grid" section above.  For Frad_SW to be correctly computed
    # on w levels, the non-constant component of F_diff and F_dir, we
    # compute taupath on w levels as a centered difference between tau
    # values on mass levels.
    #
    # This uses CLUBB's bottom-up grid ordering (k=1 at surface):
    #
    # --------taupath-->--Frad_SW----------------    k = nzt+1  (top m level)
    #        /                   \
    # -taude----------------------               k = nzt    (top t level)
    #        \                   /
    # --------taupath-->--Frad_SW----------------    k = nzt    (m level)
    #        /                   \
    # -taude----------------------               k = nzt-1  (t level)
    #        \                   /
    # --------taupath-->--Frad_SW----------------    k = nzt-1  (m level)
    #                        ...
    # --------taupath-->--Frad_SW----------------    k = 1      (sfc m level)
    #
    # At momentum level k, the thermo level above is k and below is k-1.
    # taupath accumulates from the top (k=nzt+1) downward to the surface
    # (k=1).
    #
    # Vince Larson changed the F variables to w levels  03 Feb 2005
    # Michael Falk changed the loop to start at k=2 and then solved
    # separately for k=1 so the array didn't go out of bounds.
    #
    # This code makes the same assumption as above that dwm=dmw.
    #
    # Comments by Michael Falk, 16 February 2005.                          c
    #                                                                      c
    # ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
    #                  compatible with the use of a stretched
    #                  (or unevenly-spaced) grid, as well as with the use
    #                  of an evenly-spaced grid.  Therefore, dwm is not
    #                  necessarily equal to dmw.  Interpolation functions
    #                  are used to compute any weighted averages, rather
    #                  than using general numbers such as (1/2), which is
    #                  compatible only with an evenly-spaced grid.
    #                  Brian Griffin; May 10, 2008.
    #
    # -----------------------------------------------------------------------
    # TODO(port-mirror): ``shortwave_column`` and its scan body are the local
    # JAX adaptation required to express the source's i-parallel, descending-k
    # ``taupath`` recurrence without Python scalar mutation.
    # taupath has a sequential dependency over k per column,
    # so parallelize over columns with the k sweep sequential inside.
    def shortwave_column(taude_i, zm_i, zt_i, c1_i, c2_i):
        # Top momentum level (k = nzt+1)
        taupath = jnp.where(
            l_center, 0.5 * taude_i[nzt - 1], 0.0
        )  # Running cumulative D-E optical depth          [-]
        F_diff_k = (-4.0 / 3.0) * Fs0 * (
            rp * (c1_i * jnp.exp(-rk * taupath) - c2_i * jnp.exp(rk * taupath))
            - beta * jnp.exp(-taupath / xi_abs)
        )
        F_dir_k = -Fs0 * xi_abs * jnp.exp(-taupath / xi_abs)
        Frad_SW_top = F_diff_k + F_dir_k

        # Interior momentum levels (k = nzt down to 2).
        # Going from momentum k+1 down to k, we cross thermo level k.
        # At momentum k, the thermo level above is k and below is k-1.
        def interior(taupath, k):
            interpolated = lin_interpolate_two_points(
                zm_i[k], zt_i[k], zt_i[k - 1], taude_i[k], taude_i[k - 1]
            )
            taupath = taupath + jnp.where(l_center, interpolated, taude_i[k - 1])
            F_diff_k = (-4.0 / 3.0) * Fs0 * (
                rp * (c1_i * jnp.exp(-rk * taupath) - c2_i * jnp.exp(rk * taupath))
                - beta * jnp.exp(-taupath / xi_abs)
            )
            F_dir_k = -Fs0 * xi_abs * jnp.exp(-taupath / xi_abs)
            return taupath, F_diff_k + F_dir_k

        taupath, Frad_SW_interior = jax.lax.scan(
            interior, taupath, jnp.arange(nzt - 1, 0, -1)
        )

        # Bottom momentum level (k = 1).  The ghost thermo level below the
        # surface has taude_ghost = taude(:,1), so both centered and
        # non-centered cases simply add taude(:,1).
        taupath = taupath + taude_i[0]
        F_diff_k = (-4.0 / 3.0) * Fs0 * (
            rp * (c1_i * jnp.exp(-rk * taupath) - c2_i * jnp.exp(rk * taupath))
            - beta * jnp.exp(-taupath / xi_abs)
        )
        F_dir_k = -Fs0 * xi_abs * jnp.exp(-taupath / xi_abs)
        Frad_SW_bottom = F_diff_k + F_dir_k

        Frad_SW = jnp.zeros((nzt + 1,), dtype=rcm.dtype)
        Frad_SW = Frad_SW.at[nzt].set(Frad_SW_top)
        Frad_SW = Frad_SW.at[jnp.arange(nzt - 1, 0, -1)].set(Frad_SW_interior)
        return Frad_SW.at[0].set(Frad_SW_bottom)

    return jax.vmap(shortwave_column)(taude, zm, zt, c1, c2)


__all__ = ["sunray_sw"]
