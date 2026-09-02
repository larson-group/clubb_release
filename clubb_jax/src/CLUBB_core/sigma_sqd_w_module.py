"""JAX port of ``src/CLUBB_core/sigma_sqd_w_module.F90``.

Porting deviations:
- The Fortran output argument ``sigma_sqd_w`` is returned as the function
  result, preserving functional/JIT-friendly semantics.
- Fortran loops over grid columns and levels are vectorized with JAX array
  operations.
- OpenACC data-region comments are omitted because JAX tracing handles device
  placement.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.Skx_module import Skx_func, compute_gamma_Skw
from clubb_jax.src.CLUBB_core.constants_clubb import (
    one,
    one_hundred,
    rt_tol,
    thl_tol,
    w_tol,
    w_tol_sqd,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.model_flags import l_gamma_Skw
from clubb_jax.src.CLUBB_core.grid_class import zm2zt2zm, zt2zm
from clubb_jax.src.CLUBB_core import Grid


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "l_predict_upwp_vpwp",
    ),
)
def compute_sigma_sqd_w(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    wp3,
    wp2,
    thlp2,
    rtp2,
    up2,
    vp2,
    wpthlp,
    wprtp,
    upwp,
    vpwp,
    clubb_params,
    l_predict_upwp_vpwp: bool,
):
    """Compute the variable sigma_sqd_w (PDF width parameter).

    The value of sigma_sqd_w is restricted in the ADG1 PDF in order to keep
    the marginal PDFs of all responder variables (variables that do not set
    the mixture fraction) valid.  The limits on sigma_sqd_w in order to keep
    the PDF of a responder variable, x, valid are:

    0 <= sigma_sqd_w <= 1 - <w'x'>^2 / ( <w'^2> * <x'^2> ).

    The overall limits on sigma_sqd_w must be applied based on the most
    restrictive case so that all Double Gaussian PDF responder variables, x,
    have realizable PDFs.  The overall limits on sigma_sqd_w are:

    0 <= sigma_sqd_w <= 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x).

    The equation used for sigma_sqd_w is:

    sigma_sqd_w = gamma_Skw_fnc
                  * ( 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x) );

    where 0 <= gamma_Skw_fnc <= 1.

    References:
      Eqn 22 in ``Equations for CLUBB''
    """

    # ---- Begin Code ----

    if nzm > 1:
        wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)
    else:
        # Skip interpolation if nzm == 1, this only occurs for testing purposes
        wp3_zm = wp3

    Skw_zm = Skx_func(
        nzm, ngrdcol, wp2, wp3_zm,
        w_tol, clubb_params,
    )

    gamma_Skw_fnc = compute_gamma_Skw(
        nzm, ngrdcol, Skw_zm, clubb_params,
        l_gamma_Skw,
    )

    #----------------------------------------------------------------
    # Compute sigma_sqd_w with new formula from Vince
    #----------------------------------------------------------------

    # Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
    # variables x that are Double Gaussian PDF responder variables.  This
    # includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
    # also calculated as part of the PDF, and they are included as well.
    # Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.
    max_corr_w_x_sqd = jnp.maximum(
        (
            wpthlp
            / (
                jnp.sqrt(wp2 * thlp2)
                + one_hundred * w_tol * thl_tol
            )
        ) ** 2,
        (
            wprtp
            / (
                jnp.sqrt(wp2 * rtp2)
                + one_hundred * w_tol * rt_tol
            )
        ) ** 2,
    )

    if l_predict_upwp_vpwp:
        max_corr_w_x_sqd = jnp.maximum(
            max_corr_w_x_sqd,
            jnp.maximum(
                (
                    upwp
                    / (
                        jnp.sqrt(up2 * wp2)
                        + one_hundred * w_tol_sqd
                    )
                ) ** 2,
                (
                    vpwp
                    / (
                        jnp.sqrt(vp2 * wp2)
                        + one_hundred * w_tol_sqd
                    )
                ) ** 2,
            ),
        )

    sigma_sqd_w_tmp = gamma_Skw_fnc * (one - jnp.minimum(max_corr_w_x_sqd, one))

    if nzm > 1:
        # Smooth in the vertical using interpolation.
        sigma_sqd_w = zm2zt2zm(
            nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp,
            zero_threshold,
        )
    else:
        # Skip interpolation if nzm == 1, this only occurs for testing purposes
        sigma_sqd_w = sigma_sqd_w_tmp

    return sigma_sqd_w
