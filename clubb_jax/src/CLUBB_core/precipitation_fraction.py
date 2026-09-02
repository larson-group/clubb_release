"""JAX port of precipitation_fraction.F90.

Description:
Sets overall precipitation fraction as well as the precipitation fraction
in each PDF component.

Porting deviations:
- Fortran output and inout arguments are returned so the routines remain
  functional and JIT-friendly.
- Fortran nested vertical/grid-column loops are represented with array
  expressions, reductions, and branch masks.
- Divisions in unselected ``jnp.where`` branches use safe denominators so
  JIT/AD does not carry NaNs from paths that the Fortran branch structure would
  not execute.
- ``precip_fraction`` updates immutable JAX stats state but does not perform the Fortran
  debug-level ``precip_frac_assert_check`` side effect inside the jitted path.
- ``precip_frac_assert_check`` is a NumPy validation helper that returns
  ``bool`` instead of mutating ``err_info`` and printing diagnostics.
"""

import numpy as np


from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

import jax.numpy as jnp
from jax import lax

from clubb_jax.src.CLUBB_core.constants_clubb import (
    cloud_frac_min,
    eps,
)
from clubb_jax.src.CLUBB_core.parameter_indices import iupsilon_precip_frac_rat


precip_frac_calc_type = 2
_MAX_HM_IP_COMP_MEAN = 0.0025
_PRECIP_FRAC_TOL_COEF = 0.1


def precip_fraction(
    gr,
    nzt,
    ngrdcol,
    hydromet_dim,
    hydromet,
    cloud_frac,
    cloud_frac_1,
    l_mix_rat_hm,
    l_frozen_hm,
    hydromet_tol,
    cloud_frac_2,
    ice_supersat_frac,
    ice_supersat_frac_1,
    ice_supersat_frac_2,
    mixt_frac,
    clubb_params,
    err_info,
    stats,
):
    """Determines (overall) precipitation fraction over the horizontal domain, as
    well as the precipitation fraction within each PDF component, at every
    vertical grid level.
    """
    del nzt, ngrdcol

    hydromet = jnp.asarray(hydromet, dtype=jnp.float64)
    cloud_frac = jnp.asarray(cloud_frac, dtype=jnp.float64)
    cloud_frac_1 = jnp.asarray(cloud_frac_1, dtype=jnp.float64)
    cloud_frac_2 = jnp.asarray(cloud_frac_2, dtype=jnp.float64)
    ice_supersat_frac = jnp.asarray(ice_supersat_frac, dtype=jnp.float64)
    ice_supersat_frac_1 = jnp.asarray(ice_supersat_frac_1, dtype=jnp.float64)
    ice_supersat_frac_2 = jnp.asarray(ice_supersat_frac_2, dtype=jnp.float64)
    mixt_frac = jnp.asarray(mixt_frac, dtype=jnp.float64)
    hydromet_tol = jnp.asarray(hydromet_tol, dtype=jnp.float64)
    l_mix_rat_hm = jnp.asarray(l_mix_rat_hm, dtype=bool)
    l_frozen_hm = jnp.asarray(l_frozen_hm, dtype=bool)

    # Set the minimum allowable precipitation fraction when hydrometeors are
    # found at a grid level.
    any_frozen_hm = jnp.any(l_frozen_hm)
    # Warm microphysics.
    warm_precip_frac_tol = _PRECIP_FRAC_TOL_COEF * jnp.max(cloud_frac, axis=1)
    # Ice microphysics included.
    frozen_precip_frac_tol = _PRECIP_FRAC_TOL_COEF * jnp.maximum(
        jnp.max(cloud_frac, axis=1),
        jnp.max(ice_supersat_frac, axis=1),
    )
    precip_frac_tol = jnp.maximum(
        jnp.where(any_frozen_hm, frozen_precip_frac_tol, warm_precip_frac_tol),
        cloud_frac_min,
    )
    precip_frac_tol_2d = precip_frac_tol[:, None]

    # !!! Find overall precipitation fraction.
    # The precipitation fraction is the greatest cloud fraction at or above a
    # vertical level.
    precip_frac_base = jnp.where(
        any_frozen_hm,
        jnp.maximum(cloud_frac, ice_supersat_frac),
        cloud_frac,
    )
    if int(gr.grid_dir_indx) > 0:
        precip_frac = lax.cummax(precip_frac_base, axis=1, reverse=True)
    else:
        precip_frac = lax.cummax(precip_frac_base, axis=1)

    # !!! Special checks for overall precipitation fraction
    has_hydromet = jnp.any(hydromet >= hydromet_tol[None, None, :], axis=2)
    precip_frac = jnp.where(
        has_hydromet,
        # In a scenario where we find any hydrometeor at this grid level, but
        # no cloud at or above this grid level, set precipitation fraction to
        # a minimum threshold value.
        jnp.maximum(precip_frac, precip_frac_tol_2d),
        # The means (overall) of every precipitating hydrometeor are all less
        # than their respective tolerance amounts.  They are all considered to
        # have values of 0.  There are not any hydrometeor species found at
        # this grid level.  There is also no cloud at or above this grid
        # level, so set precipitation fraction to 0.
        0.0,
    )

    # !!! Find precipitation fraction within each PDF component.
    #
    # The overall precipitation fraction, f_p, is given by the equation:
    #
    # f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    #
    # where "a" is the mixture fraction (weight of PDF component 1), f_p(1) is
    # the precipitation fraction within PDF component 1, and f_p(2) is the
    # precipitation fraction within PDF component 2.  Overall precipitation
    # fraction is found according the method above, and mixture fraction is
    # already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    # values for f_p(1) and f_p(2) must satisfy the above equation.
    if precip_frac_calc_type == 1:
        # Calculatate precip_frac_1 and precip_frac_2 based on the greatest
        # weighted cloud_frac_1 at or above a grid level.
        precip_frac_1, precip_frac_2 = component_precip_frac_weighted(
            gr,
            hydromet_dim,
            l_frozen_hm,
            hydromet_tol,
            hydromet,
            precip_frac,
            cloud_frac_1,
            cloud_frac_2,
            ice_supersat_frac_1,
            ice_supersat_frac_2,
            mixt_frac,
            precip_frac_tol,
        )
    elif precip_frac_calc_type == 2:
        # Specified method.
        clubb_params = jnp.asarray(clubb_params, dtype=jnp.float64)
        if clubb_params.ndim == 1:
            upsilon_precip_frac_rat = clubb_params[iupsilon_precip_frac_rat]
        else:
            upsilon_precip_frac_rat = clubb_params[:, iupsilon_precip_frac_rat][:, None]
        precip_frac_1, precip_frac_2 = component_precip_frac_specify(
            hydromet_dim,
            hydromet_tol,
            upsilon_precip_frac_rat,
            hydromet,
            precip_frac,
            mixt_frac,
            precip_frac_tol,
        )
    else:
        raise ValueError("Invalid option to calculate precip_frac_1 and precip_frac_2.")

    # Increase Precipiation Fraction under special conditions.
    #
    # There are scenarios that sometimes occur that require precipitation
    # fraction to be boosted.  Precipitation fraction is calculated from cloud
    # fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
    # PDF may become entirely subsaturated with respect to liquid and ice,
    # resulting in both a cloud fraction of 0 and an ice supersaturation
    # fraction of 0.  When this happens, precipitation fraction drops to 0 when
    # there aren't any hydrometeors present at that grid level, or to
    # precip_frac_tol when there is at least one hydrometeor present at that
    # grid level.  However, sometimes there are large values of hydrometeors
    # found at that grid level.  When this occurs, the PDF component in-precip
    # mean of a hydrometeor can become ridiculously large.  This is because the
    # ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
    # the equation:
    #
    # mu_hm_i = hm_i / precip_frac_i;
    #
    # where hm_i is the overall ith PDF component mean of the hydrometeor, and
    # precip_frac_i is the ith PDF component precipitation fraction.  When
    # precip_frac_i has a value of precip_frac_tol and hm_i is large, mu_hm_i
    # can be huge.  This can cause enormous microphysical process rates and
    # result in numerical instability.  It is also very inaccurate.
    #
    # In order to limit this problem, the ith PDF component precipitation
    # fraction is increased in order to decrease mu_hm_i.  First, an "upper
    # limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
    # called max_hm_ip_comp_mean.  At every vertical level and for every
    # hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
    # exceeding the "upper limit".  The check is:
    #
    # hm_i / precip_frac_i ( which = mu_hm_i )  >  max_hm_ip_comp_mean,
    #
    # which can be rewritten:
    #
    # hm_i > precip_frac_i * max_hm_ip_comp_mean.
    #
    # Since hm_i has not been calculated yet, the assumption for this check is
    # that all of the hydrometeor is found in one PDF component, which is the
    # worst-case scenario in violating this limit.  The check becomes:
    #
    # <hm> / ( mixt_frac * precip_frac_1 )  >  max_hm_ip_comp_mean;
    #    in PDF comp. 1; and
    # <hm> / ( ( 1 - mixt_frac ) * precip_frac_2 )  >  max_hm_ip_comp_mean;
    #    in PDF comp. 2.
    #
    # These limits can be rewritten as:
    #
    # <hm>  >  mixt_frac * precip_frac_1 * max_hm_ip_comp_mean;
    #    in PDF comp. 1; and
    # <hm>  >  ( 1 - mixt_frac ) * precip_frac_2 * max_hm_ip_comp_mean;
    #    in PDF comp. 2.
    #
    # When component precipitation fraction is found to be in excess of the
    # limit, precip_frac_i is increased to:
    #
    # <hm> / ( mixt_frac * max_hm_ip_comp_mean );
    #    when the limit is exceeded in PDF comp. 1; and
    # <hm> / ( ( 1 - mixt_frac ) * max_hm_ip_comp_mean );
    #    when the limit is exceeded in PDF comp. 2.
    #
    # Of course, precip_frac_i is not allowed to exceed 1, so when
    # <hm> / mixt_frac (or <hm> / ( 1 - mixt_frac )) is already greater than
    # max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
    # max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
    # compared to what it would have been using precip_frac_tol.  In the event
    # that multiple hydrometeor mixing ratios violate the check, the code is set
    # up so that precip_frac_i is increased based on the highest hm_i.
    one_minus_mixt_frac = 1.0 - mixt_frac
    mixt_frac_safe = jnp.where(mixt_frac != 0.0, mixt_frac, 1.0)
    one_minus_mixt_frac_safe = jnp.where(
        one_minus_mixt_frac != 0.0,
        one_minus_mixt_frac,
        1.0,
    )

    for ivar in range(int(hydromet_dim)):
        hydromet_i = hydromet[:, :, ivar]
        hydromet_present = hydromet_i >= hydromet_tol[ivar]
        l_mix_rat_hm_i = l_mix_rat_hm[ivar]

        # The hydrometeor is a mixing ratio.
        boost_component_1 = (
            l_mix_rat_hm_i
            & hydromet_present
            & (
                hydromet_i
                > mixt_frac * precip_frac_1 * _MAX_HM_IP_COMP_MEAN
            )
        )
        precip_frac_1_limited = jnp.maximum(
            jnp.minimum(
                hydromet_i / (mixt_frac_safe * _MAX_HM_IP_COMP_MEAN),
                1.0,
            ),
            precip_frac_tol_2d,
        )
        precip_frac_1 = jnp.where(
            boost_component_1,
            # Increase precipitation fraction in the 1st PDF component.
            precip_frac_1_limited,
            precip_frac_1,
        )

        boost_component_2 = (
            l_mix_rat_hm_i
            & hydromet_present
            & (
                hydromet_i
                > one_minus_mixt_frac
                * precip_frac_2
                * _MAX_HM_IP_COMP_MEAN
            )
        )
        precip_frac_2_limited = jnp.maximum(
            jnp.minimum(
                hydromet_i / (one_minus_mixt_frac_safe * _MAX_HM_IP_COMP_MEAN),
                1.0,
            ),
            precip_frac_tol_2d,
        )
        precip_frac_2 = jnp.where(
            boost_component_2,
            # Increase precipitation fraction in the 2nd PDF component.
            precip_frac_2_limited,
            precip_frac_2,
        )

    # Recalculate overall precipitation fraction for consistency.
    precip_frac = (
        mixt_frac * precip_frac_1
        + (1.0 - mixt_frac) * precip_frac_2
    )
    # Double check that precip_frac_tol <= precip_frac <= 1 when hydrometeors
    # are found at a grid level.
    # PLEASE DO NOT ALTER precip_frac, precip_frac_1, or precip_frac_2 anymore
    # after this point in the code.
    precip_frac = jnp.where(
        has_hydromet,
        jnp.minimum(jnp.maximum(precip_frac, precip_frac_tol_2d), 1.0),
        precip_frac,
    )

    # Statistics
    stats = stats.update("precip_frac_tol", precip_frac_tol)

    return err_info, precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol, stats


def component_precip_frac_weighted(
    gr,
    hydromet_dim,
    l_frozen_hm,
    hydromet_tol,
    hydromet,
    precip_frac,
    cloud_frac_1,
    cloud_frac_2,
    ice_supersat_frac_1,
    ice_supersat_frac_2,
    mixt_frac,
    precip_frac_tol,
):
    """Set precipitation fraction in each component of the PDF.  The weighted 1st
    PDF component precipitation fraction (weighted_pfrac_1) at a grid level is
    calculated by the greatest value of mixt_frac * cloud_frac_1 at or above
    the relevant grid level.  Likewise, the weighted 2nd PDF component
    precipitation fraction (weighted_pfrac_2) at a grid level is calculated by
    the greatest value of ( 1 - mixt_frac ) * cloud_frac_2 at or above the
    relevant grid level.

    The fraction weighted_pfrac_1 / ( weighted_pfrac_1 + weighted_pfrac_2 ) is
    the weighted_pfrac_1 fraction.  Multiplying this fraction by overall
    precipitation fraction and then dividing by mixt_frac produces the 1st PDF
    component precipitation fraction (precip_frac_1).  Then, calculate the 2nd
    PDF component precipitation fraction (precip_frac_2) accordingly.
    """
    del hydromet_dim

    hydromet = jnp.asarray(hydromet, dtype=jnp.float64)
    hydromet_tol = jnp.asarray(hydromet_tol, dtype=jnp.float64)
    precip_frac = jnp.asarray(precip_frac, dtype=jnp.float64)
    cloud_frac_1 = jnp.asarray(cloud_frac_1, dtype=jnp.float64)
    cloud_frac_2 = jnp.asarray(cloud_frac_2, dtype=jnp.float64)
    ice_supersat_frac_1 = jnp.asarray(ice_supersat_frac_1, dtype=jnp.float64)
    ice_supersat_frac_2 = jnp.asarray(ice_supersat_frac_2, dtype=jnp.float64)
    mixt_frac = jnp.asarray(mixt_frac, dtype=jnp.float64)
    precip_frac_tol = jnp.asarray(precip_frac_tol, dtype=jnp.float64)

    any_frozen_hm = jnp.any(jnp.asarray(l_frozen_hm, dtype=bool))
    one_minus_mixt_frac = 1.0 - mixt_frac
    # !!! Find precipitation fraction within PDF component 1.
    # The method used to find overall precipitation fraction will also be to
    # find precipitation fraction within PDF component 1 and PDF component 2.
    # In order to do so, it is assumed (poorly) that PDF component 1 overlaps
    # PDF component 1 and that PDF component 2 overlaps PDF component 2 at every
    # vertical level in the vertical profile.
    #
    # The weighted precipitation fraction in PDF component 1 is the greatest
    # value of the product of mixture fraction (mixt_frac) and 1st PDF
    # component cloud fraction at or above a vertical level.  Likewise, the
    # weighted precipitation fraction in PDF component 2 is the greatest
    # value of the product of ( 1 - mixt_frac ) and 2nd PDF component cloud
    # fraction at or above a vertical level.
    weighted_pfrac_1_base = jnp.where(
        any_frozen_hm,
        # Ice microphysics included.
        jnp.maximum(
            mixt_frac * cloud_frac_1,
            mixt_frac * ice_supersat_frac_1,
        ),
        # Warm microphysics.
        mixt_frac * cloud_frac_1,
    )
    weighted_pfrac_2_base = jnp.where(
        any_frozen_hm,
        # Ice microphysics included.
        jnp.maximum(
            one_minus_mixt_frac * cloud_frac_2,
            one_minus_mixt_frac * ice_supersat_frac_2,
        ),
        # Warm microphysics.
        one_minus_mixt_frac * cloud_frac_2,
    )
    if int(gr.grid_dir_indx) > 0:
        weighted_pfrac_1 = lax.cummax(weighted_pfrac_1_base, axis=1, reverse=True)
        weighted_pfrac_2 = lax.cummax(weighted_pfrac_2_base, axis=1, reverse=True)
    else:
        weighted_pfrac_1 = lax.cummax(weighted_pfrac_1_base, axis=1)
        weighted_pfrac_2 = lax.cummax(weighted_pfrac_2_base, axis=1)

    has_hydromet = jnp.any(hydromet >= hydromet_tol[None, None, :], axis=2)
    precip_frac_tol_2d = precip_frac_tol[:, None]
    mixt_frac_safe = jnp.where(mixt_frac != 0.0, mixt_frac, 1.0)
    one_minus_mixt_frac_safe = jnp.where(
        one_minus_mixt_frac != 0.0,
        one_minus_mixt_frac,
        1.0,
    )
    weighted_pfrac_sum = weighted_pfrac_1 + weighted_pfrac_2
    weighted_pfrac_sum_safe = jnp.where(
        weighted_pfrac_sum > 0.0,
        weighted_pfrac_sum,
        1.0,
    )

    # Calculate precip_frac_1 and special cases for precip_frac_1.
    # Calculate precipitation fraction in the 1st PDF component.
    precip_frac_1 = jnp.where(
        weighted_pfrac_sum > 0.0,
        # Adjust weighted 1st PDF component precipitation fraction by
        # multiplying it by a factor.  That factor is overall precipitation
        # fraction divided by the sum of the weighted 1st PDF component
        # precipitation fraction and the weighted 2nd PDF component
        # precipitation fraction.  The 1st PDF component precipitation
        # fraction is then found by dividing the adjusted weighted 1st PDF
        # component precipitation fraction by mixture fraction.
        weighted_pfrac_1
        * (precip_frac / weighted_pfrac_sum_safe)
        / mixt_frac_safe,
        # Usually, the sum of the weighted 1st PDF component precipitation
        # fraction and the 2nd PDF component precipitation fraction go to 0
        # when overall precipitation fraction goes to 0.  Since 1st PDF
        # component weighted precipitation fraction is 0, 1st PDF component
        # precipitation fraction also 0.
        0.0,
    )

    # Special cases for precip_frac_1.
    precip_frac_1_limit = jnp.minimum(1.0, precip_frac / mixt_frac_safe)
    precip_frac_1 = jnp.where(
        has_hydromet & (precip_frac_1 > precip_frac_1_limit),
        # Using the above method, it is possible for precip_frac_1 to be
        # greater than 1.  For example, the mixture fraction at level k+1 is
        # 0.10 and the cloud_frac_1 at level k+1 is 1, resulting in a
        # weighted_pfrac_1 of 0.10.  This product is greater than the product
        # of mixt_frac and cloud_frac_1 at level k.  The mixture fraction at
        # level k is 0.05, resulting in a precip_frac_1 of 2.  The value of
        # precip_frac_1 is limited at 1.  The leftover precipitation fraction
        # (a result of the decreasing weight of PDF component 1 between the
        # levels) is applied to PDF component 2.
        # Additionally, when weighted_pfrac_1 at level k is greater than
        # overall precipitation fraction at level k, the resulting calculation
        # of precip_frac_2 at level k will be negative.
        precip_frac_1_limit,
        jnp.where(
            has_hydromet
            & (precip_frac_1 > 0.0)
            & (precip_frac_1 < precip_frac_tol_2d),
            # In a scenario where we find precipitation in the 1st PDF component
            # (it is allowed to have a value of 0 when all precipitation is found
            # in the 2nd PDF component) but it is tiny (less than tolerance
            # level), boost 1st PDF component precipitation fraction to tolerance
            # level.
            precip_frac_tol_2d,
            # The means (overall) of every precipitating hydrometeor are all less
            # than their respective tolerance amounts.  They are all considered to
            # have values of 0.  There are not any hydrometeor species found at
            # this grid level.  There is also no cloud at or above this grid
            # level, so set 1st component precipitation fraction to 0.
            jnp.where(has_hydromet, precip_frac_1, 0.0),
        ),
    )

    # !!! Find precipitation fraction within PDF component 2.
    # The equation for precipitation fraction within PDF component 2 is:
    #
    # f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
    #
    # given the overall precipitation fraction, f_p (calculated above), the
    # precipitation fraction within PDF component 1, f_p(1) (calculated above),
    # and mixture fraction, a.  Any leftover precipitation fraction from
    # precip_frac_1 will be included in this calculation of precip_frac_2.
    # Calculate precipitation fraction in the 2nd PDF component.
    precip_frac_2 = jnp.where(
        has_hydromet,
        jnp.maximum(
            (precip_frac - mixt_frac * precip_frac_1)
            / one_minus_mixt_frac_safe,
            0.0,
        ),
        0.0,
    )

    # Special cases for precip_frac_2.
    precip_frac_1_if_2_gt_1 = (
        precip_frac - one_minus_mixt_frac
    ) / mixt_frac_safe
    precip_frac_2_if_1_gt_1 = (
        precip_frac - mixt_frac
    ) / one_minus_mixt_frac_safe
    precip_frac_2_if_1_lt_tol = precip_frac_tol_2d * (
        (precip_frac / precip_frac_tol_2d - mixt_frac)
        / one_minus_mixt_frac_safe
    )
    precip_frac_1_after_2_gt_1 = jnp.where(
        precip_frac_1_if_2_gt_1 > 1.0,
        1.0,
        jnp.where(
            (precip_frac_1_if_2_gt_1 > 0.0)
            & (precip_frac_1_if_2_gt_1 < precip_frac_tol_2d),
            precip_frac_tol_2d,
            precip_frac_1_if_2_gt_1,
        ),
    )
    precip_frac_2_after_2_gt_1 = jnp.where(
        precip_frac_1_if_2_gt_1 > 1.0,
        precip_frac_2_if_1_gt_1,
        jnp.where(
            (precip_frac_1_if_2_gt_1 > 0.0)
            & (precip_frac_1_if_2_gt_1 < precip_frac_tol_2d),
            precip_frac_2_if_1_lt_tol,
            1.0,
        ),
    )

    precip_frac_1_if_2_lt_tol = (
        precip_frac - one_minus_mixt_frac * precip_frac_tol_2d
    ) / mixt_frac_safe
    precip_frac_1_after_2_lt_tol = jnp.where(
        precip_frac_1_if_2_lt_tol > 1.0,
        1.0,
        jnp.where(
            (precip_frac_1_if_2_lt_tol > 0.0)
            & (precip_frac_1_if_2_lt_tol < precip_frac_tol_2d),
            precip_frac_tol_2d,
            precip_frac_1_if_2_lt_tol,
        ),
    )
    precip_frac_2_after_2_lt_tol = jnp.where(
        precip_frac_1_if_2_lt_tol > 1.0,
        precip_frac_2_if_1_gt_1,
        jnp.where(
            (precip_frac_1_if_2_lt_tol > 0.0)
            & (precip_frac_1_if_2_lt_tol < precip_frac_tol_2d),
            precip_frac_2_if_1_lt_tol,
            precip_frac_tol_2d,
        ),
    )

    precip_frac_2_gt_1 = has_hydromet & (precip_frac_2 > 1.0)
    precip_frac_2_lt_tol = (
        has_hydromet
        & (~precip_frac_2_gt_1)
        & (precip_frac_2 > 0.0)
        & (precip_frac_2 < precip_frac_tol_2d)
    )
    precip_frac_1 = jnp.where(
        precip_frac_2_gt_1,
        # Again, it is possible for precip_frac_2 to be greater than 1.
        # For example, the mixture fraction at level k+1 is 0.10 and the
        # cloud_frac_1 at level k+1 is 1, resulting in a weighted_pfrac_1
        # of 0.10.  This product is greater than the product of mixt_frac
        # and cloud_frac_1 at level k.  Additionally, precip_frac (overall)
        # is 1 for level k.  The mixture fraction at level k is 0.5,
        # resulting in a precip_frac_1 of 0.2.  Using the above equation,
        # precip_frac_2 is calculated to be 1.8.  The value of
        # precip_frac_2 is limited at 1.  The leftover precipitation
        # fraction (as a result of the increasing weight of component 1
        # between the levels) is applied to PDF component 1.
        precip_frac_1_after_2_gt_1,
        jnp.where(precip_frac_2_lt_tol, precip_frac_1_after_2_lt_tol, precip_frac_1),
    )
    precip_frac_2 = jnp.where(
        precip_frac_2_gt_1,
        precip_frac_2_after_2_gt_1,
        # In a scenario where we find precipitation in the 2nd PDF
        # component (it is allowed to have a value of 0 when all
        # precipitation is found in the 1st PDF component) but it is tiny
        # (less than tolerance level), boost 2nd PDF component
        # precipitation fraction to tolerance level.
        jnp.where(precip_frac_2_lt_tol, precip_frac_2_after_2_lt_tol, precip_frac_2),
    )

    return precip_frac_1, precip_frac_2


def component_precip_frac_specify(
    hydromet_dim,
    hydromet_tol,
    upsilon_precip_frac_rat,
    hydromet,
    precip_frac,
    mixt_frac,
    precip_frac_tol,
):
    """Calculates the precipitation fraction in each PDF component.

    The equation for precipitation fraction is:

    f_p = mixt_frac * f_p(1) + ( 1 - mixt_frac ) * f_p(2);

    where f_p is overall precipitation fraction, f_p(1) is precipitation
    fraction in the 1st PDF component, f_p(2) is precipitation fraction in the
    2nd PDF component, and mixt_frac is the mixture fraction.  Using this
    method, a new specified parameter is introduced, upsilon, where:

    upsilon = mixt_frac * f_p(1) / f_p; and where 0 <= upsilon <= 1.

    In other words, upsilon is the ratio of mixt_frac * f_p(1) to f_p.  Since
    f_p and mixt_frac are calculated previously, and upsilon is specified,
    f_p(1) can be calculated by:

    f_p(1) = upsilon * f_p / mixt_frac;

    and has an upper limit of 1.  The value of f_p(2) can then be calculated
    by:

    f_p(2) = ( f_p - mixt_frac * f_p(1) ) / ( 1 - mixt_frac );

    and also has an upper limit of 1.  When upsilon = 1, all of the
    precipitation is found in the 1st PDF component (as long as
    f_p <= mixt_frac, otherwise it would cause f_p(1) to be greater than 1).
    When upsilon = 0, all of the precipitation is found in the 2nd PDF
    component (as long as f_p <= 1 - mixt_frac, otherwise it would cause
    f_p(2) to be greater than 1).  When upsilon is between 0 and 1,
    precipitation is split between the two PDF components accordingly.
    """
    del hydromet_dim

    hydromet = jnp.asarray(hydromet, dtype=jnp.float64)
    hydromet_tol = jnp.asarray(hydromet_tol, dtype=jnp.float64)
    precip_frac = jnp.asarray(precip_frac, dtype=jnp.float64)
    mixt_frac = jnp.asarray(mixt_frac, dtype=jnp.float64)
    precip_frac_tol = jnp.asarray(precip_frac_tol, dtype=jnp.float64)[:, None]
    upsilon_precip_frac_rat = jnp.asarray(
        upsilon_precip_frac_rat,
        dtype=jnp.float64,
    )

    has_hydromet = jnp.any(hydromet >= hydromet_tol[None, None, :], axis=2)
    one_minus_mixt_frac = 1.0 - mixt_frac
    mixt_frac_safe = jnp.where(mixt_frac != 0.0, mixt_frac, 1.0)
    one_minus_mixt_frac_safe = jnp.where(
        one_minus_mixt_frac != 0.0,
        one_minus_mixt_frac,
        1.0,
    )

    # There are hydrometeors found at this grid level.
    # Loop over all vertical levels.
    precip_frac_1_one = jnp.where(
        precip_frac <= mixt_frac,
        # All the precipitation is found in the 1st PDF component.
        precip_frac / mixt_frac_safe,
        # Some precipitation is found in the 2nd PDF component.
        1.0,
    )
    precip_frac_2_one_initial = (
        precip_frac - mixt_frac
    ) / one_minus_mixt_frac_safe
    precip_frac_1_one_recalc = (
        precip_frac - one_minus_mixt_frac * precip_frac_tol
    ) / mixt_frac_safe
    precip_frac_2_one_if_1_gt_1 = (
        precip_frac - mixt_frac
    ) / one_minus_mixt_frac_safe
    precip_frac_2_one_if_1_lt_tol = precip_frac_tol * (
        (precip_frac / precip_frac_tol - mixt_frac)
        / one_minus_mixt_frac_safe
    )
    precip_frac_1_one_checked = jnp.where(
        precip_frac_1_one_recalc > 1.0,
        1.0,
        jnp.where(
            precip_frac_1_one_recalc < precip_frac_tol,
            precip_frac_tol,
            precip_frac_1_one_recalc,
        ),
    )
    precip_frac_2_one_checked = jnp.where(
        precip_frac_1_one_recalc > 1.0,
        precip_frac_2_one_if_1_gt_1,
        jnp.where(
            precip_frac_1_one_recalc < precip_frac_tol,
            # fp = a*fp1+(1-a)*fp2 solving for fp2
            precip_frac_2_one_if_1_lt_tol,
            precip_frac_tol,
        ),
    )
    precip_frac_2_one = jnp.where(
        precip_frac <= mixt_frac,
        0.0,
        jnp.where(
            (precip_frac_2_one_initial > 1.0)
            & (
                jnp.abs(precip_frac - 1.0)
                < jnp.abs(precip_frac + 1.0) / 2.0 * eps
            ),
            # Set precip_frac_2 = 1.
            1.0,
            jnp.where(
                precip_frac_2_one_initial < precip_frac_tol,
                # Since precipitation is found in the 2nd PDF component, it
                # must have a value of at least precip_frac_tol.
                precip_frac_2_one_checked,
                precip_frac_2_one_initial,
            ),
        ),
    )
    precip_frac_1_one = jnp.where(
        precip_frac <= mixt_frac,
        precip_frac_1_one,
        jnp.where(
            (precip_frac_2_one_initial > 1.0)
            & (
                jnp.abs(precip_frac - 1.0)
                < jnp.abs(precip_frac + 1.0) / 2.0 * eps
            ),
            1.0,
            jnp.where(
                precip_frac_2_one_initial < precip_frac_tol,
                # Recalculate precip_frac_1
                precip_frac_1_one_checked,
                precip_frac_1_one,
            ),
        ),
    )

    precip_frac_1_zero = jnp.where(
        precip_frac <= one_minus_mixt_frac,
        # All the precipitation is found in the 2nd PDF component.
        0.0,
        # Some precipitation is found in the 1st PDF component.
        (precip_frac - one_minus_mixt_frac) / mixt_frac_safe,
    )
    precip_frac_2_zero = jnp.where(
        precip_frac <= one_minus_mixt_frac,
        precip_frac / one_minus_mixt_frac_safe,
        1.0,
    )
    precip_frac_2_zero_recalc = (
        precip_frac - mixt_frac * precip_frac_tol
    ) / one_minus_mixt_frac_safe
    precip_frac_1_zero_if_2_gt_1 = (
        (precip_frac - 1.0) + mixt_frac
    ) / mixt_frac_safe
    precip_frac_1_zero_if_2_lt_tol = (
        (precip_frac - precip_frac_tol) / mixt_frac_safe
        + precip_frac_tol
    )
    precip_frac_2_zero_checked = jnp.where(
        precip_frac_2_zero_recalc > 1.0,
        1.0,
        jnp.where(
            precip_frac_2_zero_recalc < precip_frac_tol,
            precip_frac_tol,
            precip_frac_2_zero_recalc,
        ),
    )
    precip_frac_1_zero_checked = jnp.where(
        precip_frac_2_zero_recalc > 1.0,
        precip_frac_1_zero_if_2_gt_1,
        jnp.where(
            precip_frac_2_zero_recalc < precip_frac_tol,
            # fp = a*fp1+(1-a)*fp2 solving for fp1
            precip_frac_1_zero_if_2_lt_tol,
            precip_frac_1_zero,
        ),
    )
    precip_frac_1_zero_initial = (
        precip_frac - one_minus_mixt_frac
    ) / mixt_frac_safe
    precip_frac_1_zero = jnp.where(
        precip_frac <= one_minus_mixt_frac,
        precip_frac_1_zero,
        jnp.where(
            (precip_frac_1_zero_initial > 1.0)
            & (
                jnp.abs(precip_frac - 1.0)
                < jnp.abs(precip_frac + 1.0) / 2.0 * eps
            ),
            # Set precip_frac_1 = 1.
            1.0,
            jnp.where(
                precip_frac_1_zero_initial < precip_frac_tol,
                # Since precipitation is found in the 1st PDF component, it
                # must have a value of at least precip_frac_tol.
                precip_frac_1_zero_checked,
                precip_frac_1_zero,
            ),
        ),
    )
    precip_frac_2_zero = jnp.where(
        precip_frac <= one_minus_mixt_frac,
        precip_frac_2_zero,
        jnp.where(
            (precip_frac_1_zero_initial > 1.0)
            & (
                jnp.abs(precip_frac - 1.0)
                < jnp.abs(precip_frac + 1.0) / 2.0 * eps
            ),
            precip_frac_2_zero,
            jnp.where(
                precip_frac_1_zero_initial < precip_frac_tol,
                # Recalculate precip_frac_2
                precip_frac_2_zero_checked,
                precip_frac_2_zero,
            ),
        ),
    )

    # Precipitation is found in both PDF components.  Each component
    # must have a precipitation fraction that is at least
    # precip_frac_tol and that does not exceed 1.
    # Calculate precipitation fraction in the 1st PDF component.
    precip_frac_1_general = (
        upsilon_precip_frac_rat * precip_frac / mixt_frac_safe
    )
    # Special cases for precip_frac_1
    precip_frac_1_general = jnp.where(
        precip_frac_1_general > 1.0,
        1.0,
        jnp.where(
            precip_frac_1_general < precip_frac_tol,
            precip_frac_tol,
            precip_frac_1_general,
        ),
    )
    # Calculate precipitation fraction in the 2nd PDF component.
    precip_frac_2_general = (
        precip_frac - mixt_frac * precip_frac_1_general
    ) / one_minus_mixt_frac_safe

    # Special case for precip_frac_2
    precip_frac_1_if_2_gt_1 = (
        precip_frac - one_minus_mixt_frac
    ) / mixt_frac_safe
    precip_frac_2_if_1_gt_1 = (
        precip_frac - mixt_frac
    ) / one_minus_mixt_frac_safe
    precip_frac_2_if_1_lt_tol = precip_frac_tol * (
        (precip_frac / precip_frac_tol - mixt_frac)
        / one_minus_mixt_frac_safe
    )
    precip_frac_1_after_2_gt_1 = jnp.where(
        precip_frac_1_if_2_gt_1 > 1.0,
        1.0,
        jnp.where(
            precip_frac_1_if_2_gt_1 < precip_frac_tol,
            precip_frac_tol,
            precip_frac_1_if_2_gt_1,
        ),
    )
    precip_frac_2_after_2_gt_1 = jnp.where(
        precip_frac_1_if_2_gt_1 > 1.0,
        precip_frac_2_if_1_gt_1,
        jnp.where(
            precip_frac_1_if_2_gt_1 < precip_frac_tol,
            # fp = a*fp1+(1-a)*fp2 solving for fp2
            precip_frac_2_if_1_lt_tol,
            1.0,
        ),
    )

    precip_frac_1_if_2_lt_tol = (
        precip_frac - one_minus_mixt_frac * precip_frac_tol
    ) / mixt_frac_safe
    precip_frac_1_after_2_lt_tol = jnp.where(
        precip_frac_1_if_2_lt_tol > 1.0,
        1.0,
        jnp.where(
            precip_frac_1_if_2_lt_tol < precip_frac_tol,
            precip_frac_tol,
            precip_frac_1_if_2_lt_tol,
        ),
    )
    precip_frac_2_after_2_lt_tol = jnp.where(
        precip_frac_1_if_2_lt_tol > 1.0,
        precip_frac_2_if_1_gt_1,
        jnp.where(
            precip_frac_1_if_2_lt_tol < precip_frac_tol,
            precip_frac_2_if_1_lt_tol,
            precip_frac_tol,
        ),
    )
    precip_frac_1_general = jnp.where(
        precip_frac_2_general > 1.0,
        # Set precip_frac_2 to 1.
        # Recalculate precipitation fraction in the 1st PDF component.
        precip_frac_1_after_2_gt_1,
        jnp.where(
            precip_frac_2_general < precip_frac_tol,
            # Set precip_frac_2 to precip_frac_tol.
            # Recalculate precipitation fraction in the 1st PDF component.
            precip_frac_1_after_2_lt_tol,
            precip_frac_1_general,
        ),
    )
    precip_frac_2_general = jnp.where(
        precip_frac_2_general > 1.0,
        precip_frac_2_after_2_gt_1,
        jnp.where(
            precip_frac_2_general < precip_frac_tol,
            precip_frac_2_after_2_lt_tol,
            precip_frac_2_general,
        ),
    )

    l_upsilon_one = (
        jnp.abs(upsilon_precip_frac_rat - 1.0)
        < jnp.abs(upsilon_precip_frac_rat + 1.0) / 2.0 * eps
    )
    l_upsilon_zero = (
        jnp.abs(upsilon_precip_frac_rat)
        < jnp.abs(upsilon_precip_frac_rat) / 2.0 * eps
    )
    precip_frac_1 = jnp.where(
        l_upsilon_one,
        precip_frac_1_one,
        jnp.where(l_upsilon_zero, precip_frac_1_zero, precip_frac_1_general),
    )
    precip_frac_2 = jnp.where(
        l_upsilon_one,
        precip_frac_2_one,
        jnp.where(l_upsilon_zero, precip_frac_2_zero, precip_frac_2_general),
    )

    # There aren't any hydrometeors found at the grid level.
    precip_frac_1 = jnp.where(has_hydromet, precip_frac_1, 0.0)
    precip_frac_2 = jnp.where(has_hydromet, precip_frac_2, 0.0)

    return precip_frac_1, precip_frac_2


def precip_frac_assert_check(
    hydromet,
    hydromet_tol,
    mixt_frac,
    precip_frac,
    precip_frac_1,
    precip_frac_2,
    precip_frac_tol,
):
    """Assertion check for the precipitation fraction code."""
    hydromet = np.asarray(hydromet, dtype=np.float64)
    hydromet_tol = np.asarray(hydromet_tol, dtype=np.float64)
    mixt_frac = np.asarray(mixt_frac, dtype=np.float64)
    precip_frac = np.asarray(precip_frac, dtype=np.float64)
    precip_frac_1 = np.asarray(precip_frac_1, dtype=np.float64)
    precip_frac_2 = np.asarray(precip_frac_2, dtype=np.float64)
    precip_frac_tol = float(np.asarray(precip_frac_tol, dtype=np.float64))

    # Loop over all vertical levels.
    has_hydromet = np.any(hydromet >= hydromet_tol, axis=-1)
    # Overall precipitation fraction cannot be less than precip_frac_tol
    # when a hydrometeor is present at a grid level.
    # Overall precipitation fraction cannot exceed 1.
    # Precipitation fraction in the 1st PDF component is allowed to be 0
    # when all the precipitation is found in the 2nd PDF component.
    # Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
    # is present at a grid level.  In other words, it cannot have a value
    # that is greater than 0 but less than precip_frac_tol
    # Precipitation fraction in the 1st PDF component cannot exceed 1.
    # Precipiation fraction in the 1st PDF component cannot be negative.
    # Precipitation fraction in the 2nd PDF component is allowed to be 0
    # when all the precipitation is found in the 1st PDF component.
    # Otherwise, it cannot be less than precip_frac_tol when a hydrometeor
    # is present at a grid level.  In other words, it cannot have a value
    # that is greater than 0 but less than precip_frac_tol
    # Precipitation fraction in the 2nd PDF component cannot exceed 1.
    # Precipiation fraction in the 2nd PDF component cannot be negative.
    cloudy_ok = (
        (precip_frac >= precip_frac_tol)
        & (precip_frac <= 1.0)
        & ~((precip_frac_1 > 0.0) & (precip_frac_1 < precip_frac_tol - eps))
        & (precip_frac_1 >= 0.0)
        & (precip_frac_1 <= 1.0)
        & ~((precip_frac_2 > 0.0) & (precip_frac_2 < precip_frac_tol - eps))
        & (precip_frac_2 >= 0.0)
        & (precip_frac_2 <= 1.0)
    )
    # Overall precipitation fraction must be 0 when no hydrometeors are
    # found at a grid level.
    # Precipitation fraction in the 1st PDF component must be 0 when no
    # hydrometeors are found at a grid level.
    # Precipitation fraction in the 2nd PDF component must be 0 when no
    # hydrometeors are found at a grid level.
    clear_ok = (
        (np.abs(precip_frac) <= eps)
        & (np.abs(precip_frac_1) <= eps)
        & (np.abs(precip_frac_2) <= eps)
    )
    # The precipitation fraction equation is:
    #
    # precip_frac
    #    = mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2;
    #
    # which means that:
    #
    # precip_frac
    # - ( mixt_frac * precip_frac_1 + ( 1 - mixt_frac ) * precip_frac_2 )
    # = 0.
    #
    # Check that this is true with numerical round off.
    mixture_ok = (
        precip_frac
        - (mixt_frac * precip_frac_1 + (1.0 - mixt_frac) * precip_frac_2)
    ) <= eps

    return bool(np.all(np.where(has_hydromet, cloudy_ok, clear_ok)) and np.all(mixture_ok))


__all__ = [
    "precip_frac_calc_type",
    "precip_fraction",
    "component_precip_frac_weighted",
    "component_precip_frac_specify",
    "precip_frac_assert_check",
]
