"""JAX implementations of selected routines from ``Nc_Ncn_eqns.F90``.

Description:
Equations are provided to perform calculations back-and-forth between Nc and
Ncn, where Nc is cloud droplet concentration and Ncn is simplified cloud
nuclei concentration.  The equation that relates the two is:

Nc = Ncn * H(chi);

where chi is extended liquid water mixing ratio, which is equal to cloud
water mixing ratio, rc, when both are positive.  However, chi is negative in
subsaturated air.

Equation are provided relating mean cloud droplet concentration (overall),
Ncm, and/or mean cloud droplet concentration (in-cloud), Nc_in_cloud, to
mean simplified cloud nuclei concentration, Ncnm.

Notes:

Meaning of Nc flag combinations:

l_const_Nc_in_cloud:
When this flag is enabled, cloud droplet concentration (in-cloud) is
constant (spatially) at a grid level (it is constant over the subgrid
domain, but could vary over time depending on the value of l_predict_Nc).
The value of in-cloud Nc does not vary at a grid level.  This also means
that Ncn is constant across the entire grid level.  When this flag is turned
off, both in-cloud Nc and Ncn vary at a grid level.

l_predict_Nc:
When this flag is enabled, Nc_in_cloud (or alternatively Ncm) is predicted.
It is advanced every time step by a predictive equation, and can change
at every time step at a grid level.  When this flag is turned off,
Nc_in_cloud does not change at a grid level over the course of a model run.

References:

Porting deviations:
Fortran elemental routines are represented as broadcasting JAX functions.
Fortran branch bodies execute conditionally; JAX precomputes values for
``where`` selection, so guarded denominators are used to avoid invalid work in
unselected branches.  The long Fortran flow-chart module comment describes
full CLUBB microphysics orchestration outside these local conversion helpers
and is not repeated here.
"""

import jax.scipy.special as jsp
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    chi_tol,
    cloud_frac_min,
    eps,
    Ncn_tol,
    sqrt_2,
)


def Ncnm_to_Nc_in_cloud(
    mu_chi_1,
    mu_chi_2,
    mu_Ncn_1,
    mu_Ncn_2,
    sigma_chi_1,
    sigma_chi_2,
    sigma_Ncn_1,
    sigma_Ncn_2,
    sigma_Ncn_1_n,
    sigma_Ncn_2_n,
    corr_chi_Ncn_1_n,
    corr_chi_Ncn_2_n,
    mixt_frac,
    cloud_frac_1,
    cloud_frac_2,
):
    """Calculate Nc_in_cloud from Ncn PDF parameters.

    Description:
    The in-cloud mean of cloud droplet concentration is calculated from the
    PDF parameters involving simplified cloud nuclei concentration, Ncn, and
    cloud fraction.  At any point, cloud droplet concentration, Nc, is given
    by:

    Nc = Ncn * H(chi);

    where extended liquid water mixing ratio, chi, is equal to cloud water
    ratio, rc, when positive.  When the atmosphere is saturated at this point,
    cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    and Nc = 0.

    The overall mean of cloud droplet concentration, <Nc>, is calculated from
    the PDF parameters involving Ncn.  The in-cloud mean of cloud droplet
    concentration is calculated from <Nc> and cloud fraction.

    References:
    """
    # Calculate overall cloud fraction as calculated by the PDF.
    # The variable cloud_frac is not used here because it is altered by factors
    # such as the trapezoidal rule calculation.
    # Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
    # as long neither of these variables are altered by any factor.  They can
    # only be calculated from PDF.
    cloud_frac = mixt_frac * cloud_frac_1 + (1.0 - mixt_frac) * cloud_frac_2

    # There is cloud found at this grid level.  Calculate Nc_in_cloud.
    Ncm = Ncnm_to_Ncm(
        mu_chi_1,
        mu_chi_2,
        mu_Ncn_1,
        mu_Ncn_2,
        sigma_chi_1,
        sigma_chi_2,
        sigma_Ncn_1,
        sigma_Ncn_2,
        sigma_Ncn_1_n,
        sigma_Ncn_2_n,
        corr_chi_Ncn_1_n,
        corr_chi_Ncn_2_n,
        mixt_frac,
    )

    cloud_frac_safe = jnp.where(cloud_frac > cloud_frac_min, cloud_frac, 1.0)

    # This level is entirely clear.  Set Nc_in_cloud to <Ncn>.
    # Since <Ncn> = mu_Ncn_1 = mu_Ncn_2, use mu_Ncn_1 here.
    return jnp.where(cloud_frac > cloud_frac_min, Ncm / cloud_frac_safe, mu_Ncn_1)


def Nc_in_cloud_to_Ncnm(
    mu_chi_1,
    mu_chi_2,
    sigma_chi_1,
    sigma_chi_2,
    mixt_frac,
    Nc_in_cloud,
    cloud_frac_1,
    cloud_frac_2,
    const_Ncnp2_on_Ncnm2,
    const_corr_chi_Ncn,
):
    """Calculate Ncnm from Nc_in_cloud and PDF parameters.

    Description:
    The overall mean of simplified cloud nuclei concentration, <Ncn>, is
    calculated from the in-cloud mean of cloud droplet concentration, <Nc>,
    cloud fraction, and some of the PDF parameters.

    At any point, cloud droplet concentration, Nc, is given by:

    Nc = Ncn * H(chi);

    where extended liquid water mixing ratio, chi, is equal to cloud water
    ratio, rc, when positive.  When the atmosphere is saturated at this point,
    cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    and Nc = 0.

    The overall mean of cloud droplet concentration, <Nc>, is calculated from
    Nc_in_cloud and cloud fraction.  The value of <Ncn> is calculated from
    <Nc> and PDF parameters.

    References:
    """
    # Calculate overall cloud fraction as calculated by the PDF.
    # The variable cloud_frac is not used here because it is altered by factors
    # such as the trapezoidal rule calculation.
    # Cloud fraction can be recalculated here from cloud_frac_1 and cloud_frac_2
    # as long neither of these variables are altered by any factor.  They can
    # only be calculated from the PDF.
    cloud_frac = mixt_frac * cloud_frac_1 + (1.0 - mixt_frac) * cloud_frac_2

    # When Ncn is constant a a grid level, it is equal to Nc_in_cloud.
    # Additionally, when a level is entirely clear, <Ncn>, which is based on
    # Nc_in_cloud, here, must be set to something.  Set <Ncn> to Nc_in_cloud.
    Ncnm = Nc_in_cloud

    # There is cloud found at this grid level.  Additionally, Ncn varies.
    # Calculate Nc_in_cloud.
    Ncm = Nc_in_cloud * cloud_frac

    Ncnm_varying = Ncm_to_Ncnm(
        mu_chi_1,
        mu_chi_2,
        sigma_chi_1,
        sigma_chi_2,
        mixt_frac,
        Ncm,
        const_Ncnp2_on_Ncnm2,
        const_corr_chi_Ncn,
        Nc_in_cloud,
    )

    l_varying_Ncn = (
        (cloud_frac > cloud_frac_min)
        & (jnp.abs(const_corr_chi_Ncn * const_Ncnp2_on_Ncnm2) > eps)
    )
    return jnp.where(l_varying_Ncn, Ncnm_varying, Ncnm)


def Ncnm_to_Ncm(
    mu_chi_1,
    mu_chi_2,
    mu_Ncn_1,
    mu_Ncn_2,
    sigma_chi_1,
    sigma_chi_2,
    sigma_Ncn_1,
    sigma_Ncn_2,
    sigma_Ncn_1_n,
    sigma_Ncn_2_n,
    corr_chi_Ncn_1_n,
    corr_chi_Ncn_2_n,
    mixt_frac,
):
    """Calculate Ncm from Ncn PDF parameters.

    Description:
    The overall mean of cloud droplet concentration, <Nc>, is calculated from
    the PDF parameters involving the simplified cloud nuclei concentration,
    Ncn.  At any point, cloud droplet concentration, Nc, is given by:

    Nc = Ncn * H(chi);

    where extended liquid water mixing ratio, chi, is equal to cloud water
    ratio, rc, when positive.  When the atmosphere is saturated at this point,
    cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    and Nc = 0.

    The overall mean of cloud droplet concentration, <Nc>, is found by
    integrating over the PDF of chi and Ncn, such that:

    <Nc> = INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P(chi,Ncn) dNcn dchi;

    which can also be written as:

    <Nc> = SUM(i=1,n) mixt_frac_i
           * INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;

    where n is the number of multivariate joint PDF components, mixt_frac_i is
    the weight of the ith PDF component, and P_i is the functional form of the
    multivariate joint PDF in the ith PDF component.

    This equation is rewritten as:

    <Nc> = SUM(i=1,n) mixt_frac_i
           * INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi.

    References:
    """
    # Calculate mean cloud droplet concentration (overall), <Nc>.
    Ncm = (
        mixt_frac
        * bivar_NL_chi_Ncn_mean(
            mu_chi_1,
            mu_Ncn_1,
            sigma_chi_1,
            sigma_Ncn_1,
            sigma_Ncn_1_n,
            corr_chi_Ncn_1_n,
        )
        + (1.0 - mixt_frac)
        * bivar_NL_chi_Ncn_mean(
            mu_chi_2,
            mu_Ncn_2,
            sigma_chi_2,
            sigma_Ncn_2,
            sigma_Ncn_2_n,
            corr_chi_Ncn_2_n,
        )
    )

    return Ncm


def Ncm_to_Ncnm(
    mu_chi_1,
    mu_chi_2,
    sigma_chi_1,
    sigma_chi_2,
    mixt_frac,
    Ncm,
    const_Ncnp2_on_Ncnm2,
    const_corr_chi_Ncn,
    Ncnm_val_denom_0,
):
    """Calculate Ncnm from Ncm and PDF parameters.

    Description:
    The overall mean of simplified cloud nuclei concentration, <Ncn>, is
    calculated from the overall mean of cloud droplet concentration, <Nc>, and
    some of the PDF parameters.

    At any point, cloud droplet concentration, Nc, is given by:

    Nc = Ncn * H(chi);

    where extended liquid water mixing ratio, chi, is equal to cloud water
    ratio, rc, when positive.  When the atmosphere is saturated at this point,
    cloud water is found, and Nc = Ncn.  Otherwise, only clear air is found,
    and Nc = 0.

    Solving for <Ncn>
    =================

    In order to isolate <Ncn>, the value of <Ncn'^2>/<Ncn>^2 is set to a
    constant value, const_Ncn.  The value of this constant does not depend on
    <Ncn>.  Likewise, the value of rho_chi_Ncn does not depend on <Ncn>.
    Solving for <Ncn>, the equation becomes:

    <Ncn>
    = <Nc> / ( SUM(i=1,n) mixt_frac_i
                 ---
                 | (1/2) * erfc( - ( 1 / sqrt(2) )
                 |                 * ( ( mu_chi_i / sigma_chi_i )
                 |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
                 | where sigma_chi_i > 0 and const_Ncn > 0;
                 |
               * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
                 | where sigma_chi_i > 0 and const_Ncn = 0;
                 |
                 | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
                 |
                 | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
                 ---

    When the denominator term is 0, there is only clear air.  Both the
    numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
    to an appropriate value.

    References:
    """
    # Denominator in the equation for <Ncn>
    denominator = (
        mixt_frac
        * bivar_Ncnm_eqn_comp(
            mu_chi_1,
            sigma_chi_1,
            const_Ncnp2_on_Ncnm2,
            const_corr_chi_Ncn,
        )
        + (1.0 - mixt_frac)
        * bivar_Ncnm_eqn_comp(
            mu_chi_2,
            sigma_chi_2,
            const_Ncnp2_on_Ncnm2,
            const_corr_chi_Ncn,
        )
    )

    denominator_safe = jnp.where(denominator > 0.0, denominator, 1.0)

    # When the denominator is 0, it is usually because there is only clear
    # air.  In that scenario, Ncm should also be 0.  Set Ncnm to a value that
    # is usual or typical
    return jnp.where(denominator > 0.0, Ncm / denominator_safe, Ncnm_val_denom_0)


def bivar_NL_chi_Ncn_mean(
    mu_chi_i,
    mu_Ncn_i,
    sigma_chi_i,
    sigma_Ncn_i,
    sigma_Ncn_i_n,
    corr_chi_Ncn_i_n,
):
    """Evaluate the per-component normal-lognormal Nc integral.

    Description:
    The double integral over Ncn * H(chi) multiplied by the
    bivariate normal-lognormal joint PDF of chi and Ncn is evaluated.  The
    integral is given by:

    INT(-inf:inf) INT(0:inf) Ncn * H(chi) * P_i(chi,Ncn) dNcn dchi;

    which reduces to:

    INT(0:inf) INT(0:inf) Ncn * P_i(chi,Ncn) dNcn dchi;

    where the individual marginal distribution of chi is normal in the ith PDF
    component and the individual marginal distribution of Ncn is lognormal in
    the ith PDF component.

    References:
    """
    chi_and_Ncn_constant = (sigma_chi_i <= chi_tol) & (sigma_Ncn_i <= Ncn_tol)
    chi_constant = sigma_chi_i <= chi_tol
    Ncn_constant = sigma_Ncn_i <= Ncn_tol

    # The ith PDF component variances of both chi and Ncn are 0.
    # The ith PDF component variance of chi is 0.
    constant_value = jnp.where(mu_chi_i > 0.0, mu_Ncn_i, 0.0)

    sigma_chi_i_safe = jnp.where(sigma_chi_i > chi_tol, sigma_chi_i, 1.0)

    # The ith PDF component variance of Ncn is 0.
    Ncn_constant_value = (
        mu_Ncn_i
        * 0.5
        * jsp.erfc(-(mu_chi_i / (sqrt_2 * sigma_chi_i_safe)))
    )

    # Both chi and Ncn vary in the ith PDF component.
    both_vary_value = (
        0.5
        * mu_Ncn_i
        * jsp.erfc(
            -(1.0 / sqrt_2)
            * ((mu_chi_i / sigma_chi_i_safe) + corr_chi_Ncn_i_n * sigma_Ncn_i_n)
        )
    )

    return jnp.where(
        chi_and_Ncn_constant,
        constant_value,
        jnp.where(
            chi_constant,
            constant_value,
            jnp.where(Ncn_constant, Ncn_constant_value, both_vary_value),
        ),
    )


def bivar_Ncnm_eqn_comp(
    mu_chi_i,
    sigma_chi_i,
    const_Ncnp2_on_Ncnm2,
    const_corr_chi_Ncn,
):
    """Calculate one PDF component's denominator term in the Ncnm equation.

    Description:
    When <Ncn> is found based on the value of <Nc>, the following equation is
    used:

    <Ncn>
    = <Nc> / ( SUM(i=1,n) mixt_frac_i
                 ---
                 | (1/2) * erfc( - ( 1 / sqrt(2) )
                 |                 * ( ( mu_chi_i / sigma_chi_i )
                 |                     + rho_chi_Ncn * sqrt( const_Ncn ) ) );
                 | where sigma_chi_i > 0 and const_Ncn > 0;
                 |
               * | (1/2) * erfc( - ( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
                 | where sigma_chi_i > 0 and const_Ncn = 0;
                 |
                 | 1; where sigma_chi_i = 0 and mu_chi_i > 0;
                 |
                 | 0; where sigma_chi_i = 0 and mu_chi_i <= 0               ).
                 ---

    In the above equation, const_Ncn = <Ncn'^2> / <Ncn>^2.  It is a constant,
    prescribed parameter.  Likewise, rho_chi_Ncn is a parameter that is not
    based on the value of <Ncn>.

    When the denominator term is 0, there is only clear air.  Both the
    numerator (<Nc>) and the denominator have a value of 0, and <Ncn> is set
    to an appropriate value.

    The contribution of the ith PDF component to the denominator term in the
    equation is calculated here.

    References:
    """
    # The ith PDF component variances of chi is 0.  The value of the ith PDF
    # component variance of Ncn does not matter in this scenario.
    constant_value = jnp.where(mu_chi_i > 0.0, 1.0, 0.0)
    sigma_chi_i_safe = jnp.where(sigma_chi_i > chi_tol, sigma_chi_i, 1.0)

    # Both chi and Ncn vary in the ith PDF component.
    varying_value = (
        0.5
        * jsp.erfc(
            -(1.0 / sqrt_2)
            * (
                (mu_chi_i / sigma_chi_i_safe)
                + const_corr_chi_Ncn * jnp.sqrt(const_Ncnp2_on_Ncnm2)
            )
        )
    )

    return jnp.where(sigma_chi_i <= chi_tol, constant_value, varying_value)
