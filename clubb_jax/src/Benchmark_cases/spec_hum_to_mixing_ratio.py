"""This module contains subroutines that take various moisture variables that
were given in terms of specific humidity and converts them to terms of mixing
ratio.
"""

from __future__ import annotations


# =============================================================================
def flux_spec_hum_to_mixing_ratio(ngrdcol, rtm_zm, wpqtp):
    """This subroutine takes a flux given in terms of specific humidity,
    w'q_t', and converts it to a flux given in terms of mixing ratio, w'r_t'.

    This relationship is given by the equation:

    w'r_t' = ( 1 + r_tm )^2 * w'q_t';

    where r_tm is the mean total water mixing ratio. Higher-order terms, which
    were comparatively small in magnitude, were neglected for purposes of
    linearizing the equation.

    This subroutine has been written very generally, in order to allow for this
    conversion at any vertical level. However, the most likely application of
    this is at the surface, taking a specified surface flux in terms of
    specific humidity, w'q_t'|_sfc, and converting it to a surface flux in
    terms of mixing ratio, w'r_t'|_sfc.
    """
    del ngrdcol

    # Solve for flux in terms of total water mixing ratio.
    return (1.0 + rtm_zm) ** 2 * wpqtp


# =============================================================================
def force_spec_hum_to_mixing_ratio(ngrdcol, nzt, rtm, qtm_forcing):
    """This subroutine takes a forcing given in terms of specific humidity,
    d(q_tm)/dt|_f, and converts it to a forcing given in terms of mixing ratio,
    d(r_tm)/dt|_f.

    This relationship is given by the equation:

    d(r_tm)/dt|_f = ( 1 + r_tm )^2 * d(q_tm)/dt|_f;

    where r_tm is the mean total water mixing ratio.
    """
    del ngrdcol, nzt

    # Solve for forcing in terms of total water mixing ratio.
    return (1.0 + rtm) ** 2 * qtm_forcing


__all__ = ["flux_spec_hum_to_mixing_ratio", "force_spec_hum_to_mixing_ratio"]
