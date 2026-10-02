"""Hydrometeor mixed moments, mirroring mixed_moment_PDF_integrals.F90.

JAX adaptation: thermodynamic-level and column loops are batched arrays;
static hydrometeor loops and source routine ordering are retained.
"""
import math
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.pdf_utilities import compute_mean_binormal, calc_corr_rt_x, calc_corr_thl_x
from clubb_jax.src.CLUBB_core.constants_clubb import rt_tol, thl_tol, w_tol
from clubb_jax.src.CLUBB_core.index_mapping import hydromet2pdf_idx
from clubb_jax.src.CLUBB_core.grid_class import zt2zm


def hydrometeor_mixed_moments(gr, ngrdcol, nzt, pdf_dim, hydromet_dim,
    hydromet, hm_metadata, mu_x_1_n, mu_x_2_n, sigma_x_1_n,
    sigma_x_2_n, corr_array_1_n, corr_array_2_n, pdf_params, hydromet_pdf_params,
    precip_fracs, stats):
    # Description:
    # Calculates <rt'hm'>, <thl'hm'>, and <w'^2 hm'>, for all hydrometeors, hm.
    # These terms are used in the liquid/ice water loading term as part of the
    # buoyancy term in some of the CLUBB predictive equations.
    # References:
    #-----------------------------------------------------------------------
    rtphmp_zt = jnp.zeros_like(hydromet)
    thlphmp_zt = jnp.zeros_like(hydromet)
    wp2hmp = jnp.zeros_like(hydromet)
    hmxphmyp_zt = jnp.zeros(hydromet.shape + (hydromet_dim,))
    # Loop over all thermodynamic levels between the model lower and upper
    # boundaries (thermodynamic levels 1 to gr%nzt).


    # Unpack the means of w, rt, and thl in each PDF component.
    mu_w_1   = mu_x_1_n[..., hm_metadata.iiPDF_w]
    mu_w_2   = mu_x_2_n[..., hm_metadata.iiPDF_w]
    mu_rt_1  = pdf_params.rt_1
    mu_rt_2  = pdf_params.rt_2
    mu_thl_1 = pdf_params.thl_1
    mu_thl_2 = pdf_params.thl_2

    # Unpack the standard deviations of w, rt, and thl in each PDF component.
    sigma_w_1   = sigma_x_1_n[..., hm_metadata.iiPDF_w]
    sigma_w_2   = sigma_x_2_n[..., hm_metadata.iiPDF_w]
    sigma_rt_1  = jnp.sqrt( pdf_params.varnce_rt_1 )
    sigma_rt_2  = jnp.sqrt( pdf_params.varnce_rt_2 )
    sigma_thl_1 = jnp.sqrt( pdf_params.varnce_thl_1 )
    sigma_thl_2 = jnp.sqrt( pdf_params.varnce_thl_2 )

    # Unpack the standard deviations of chi and eta in each PDF component.
    sigma_chi_1 = sigma_x_1_n[..., hm_metadata.iiPDF_chi]
    sigma_chi_2 = sigma_x_2_n[..., hm_metadata.iiPDF_chi]
    sigma_eta_1 = sigma_x_1_n[..., hm_metadata.iiPDF_eta]
    sigma_eta_2 = sigma_x_2_n[..., hm_metadata.iiPDF_eta]

    # Unpack the mixture fraction.
    mixt_frac = pdf_params.mixt_frac

    # Unpack the precipitation fraction in each PDF component.
    precip_frac_1 = precip_fracs.precip_frac_1
    precip_frac_2 = precip_fracs.precip_frac_2

    # Unpack the coefficients of rt and thl in the chi/eta PDF transformation
    # equations for each PDF component.
    crt_1  = pdf_params.crt_1
    crt_2  = pdf_params.crt_2
    cthl_1 = pdf_params.cthl_1
    cthl_2 = pdf_params.cthl_2

    # Re-calculate rtm, thlm, and wm from PDF parameters.
    # This needs to be done because rtm and thlm have been advanced since
    # the PDF parameters have been calculated.  It is necessary to use values
    # of the mean fields (rtm, thlm, and wm) that are consistent with the
    # PDF parameters.  This does not need to be done for hydromet because
    # hydrometeors have not been advanced since the hydrometeor PDF
    # parameters were set up.
    wm   = compute_mean_binormal( mu_w_1, mu_w_2, mixt_frac )
    rtm  = compute_mean_binormal( mu_rt_1, mu_rt_2, mixt_frac )
    thlm = compute_mean_binormal( mu_thl_1, mu_thl_2, mixt_frac )


    # Calculate <rt'hm'>, <thl'hm'>, and <w'^2 hm'> for each hydrometeor
    # species.
    for hm_idx in range(hydromet_dim):

        # Unpack the mean (in-precip) of hm in each PDF component.
        mu_hm_1 = hydromet_pdf_params.mu_hm_1[..., hm_idx]
        mu_hm_2 = hydromet_pdf_params.mu_hm_2[..., hm_idx]

        # Unpack the standard deviation (in-precip) of hm in each PDF
        # component.
        sigma_hm_1 = hydromet_pdf_params.sigma_hm_1[..., hm_idx]
        sigma_hm_2 = hydromet_pdf_params.sigma_hm_2[..., hm_idx]

        # Calculate the correlation (in-precip) of rt/thl and hm for each PDF
        # component.  Since CLUBB uses a PDF transformation from rt and
        # theta-l coordinates to chi and eta coordinates for each PDF
        # component, the correlation arrays are written in terms of chi and
        # eta correlations.  This makes a calculation necessary for these
        # correlations.
        corr_chi_hm_1 = hydromet_pdf_params.corr_chi_hm_1[..., hm_idx]
        corr_chi_hm_2 = hydromet_pdf_params.corr_chi_hm_2[..., hm_idx]
        corr_eta_hm_1 = hydromet_pdf_params.corr_eta_hm_1[..., hm_idx]
        corr_eta_hm_2 = hydromet_pdf_params.corr_eta_hm_2[..., hm_idx]

        corr_rt_hm_1 = calc_corr_rt_x( crt_1, sigma_rt_1, sigma_chi_1, sigma_eta_1, corr_chi_hm_1, corr_eta_hm_1 )

        corr_rt_hm_2 = calc_corr_rt_x( crt_2, sigma_rt_2, sigma_chi_2, sigma_eta_2, corr_chi_hm_2, corr_eta_hm_2 )

        corr_thl_hm_1 = calc_corr_thl_x( cthl_1, sigma_thl_1, sigma_chi_1, sigma_eta_1, corr_chi_hm_1, corr_eta_hm_1 )

        corr_thl_hm_2 = calc_corr_thl_x( cthl_2, sigma_thl_2, sigma_chi_2, sigma_eta_2, corr_chi_hm_2, corr_eta_hm_2 )

        # Unpack the tolerance value for the hydrometeor, hm.
        hm_tol = hm_metadata.hydromet_tol[hm_idx]

        # Calculate <rt'hm'>.
        rtphmp_zt = rtphmp_zt.at[..., hm_idx].set(xphmp_integral_covar(
            mu_rt_1, mu_rt_2, mu_hm_1, mu_hm_2,
            sigma_rt_1, sigma_rt_2, sigma_hm_1, sigma_hm_2,
            corr_rt_hm_1, corr_rt_hm_2, mixt_frac, precip_frac_1,
            precip_frac_2, rtm, rt_tol, hm_tol))

        # Calculate <thl'hm'>.
        thlphmp_zt = thlphmp_zt.at[..., hm_idx].set(xphmp_integral_covar(
            mu_thl_1, mu_thl_2, mu_hm_1, mu_hm_2,
            sigma_thl_1, sigma_thl_2, sigma_hm_1, sigma_hm_2,
            corr_thl_hm_1, corr_thl_hm_2, mixt_frac, precip_frac_1,
            precip_frac_2, thlm, thl_tol, hm_tol))

        # Find the index of hydrometeor in the PDF indices.
        pdf_idx = hydromet2pdf_idx(hm_idx,hm_metadata)

        # Unpack the mean (in-precip) of ln hm in each PDF component.
        mu_hm_1_n = mu_x_1_n[..., pdf_idx]
        mu_hm_2_n = mu_x_2_n[..., pdf_idx]

        # Unpack the standard deviation (in-precip) of ln hm in each PDF
        # component.
        sigma_hm_1_n = sigma_x_1_n[..., pdf_idx]
        sigma_hm_2_n = sigma_x_2_n[..., pdf_idx]

        # Unpack the correlation (in-precip) of w and ln hm in each PDF
        # component.
        corr_w_hm_1_n = corr_array_1_n[..., pdf_idx,hm_metadata.iiPDF_w]
        corr_w_hm_2_n = corr_array_2_n[..., pdf_idx,hm_metadata.iiPDF_w]

        # Unpack the mean (overall) value of the hydrometeor.
        hm_mean = hydromet[..., hm_idx]

        # The general form of the mixed moment equation is <w'^a hm'^b>.
        # For <w'^2 hm'>, a = 2 and b = 1.
        a_exp = 2
        b_exp = 1

        # Calculate <w'^2 hm'>.
        wp2hmp = wp2hmp.at[..., hm_idx].set(xp_a_hmpb_integrals_all_MM(
            mu_w_1, mu_w_2, mu_hm_1, mu_hm_2, mu_hm_1_n, mu_hm_2_n,
            sigma_w_1, sigma_w_2, sigma_hm_1, sigma_hm_2, sigma_hm_1_n, sigma_hm_2_n,
            corr_w_hm_1_n, corr_w_hm_2_n, mixt_frac, precip_frac_1, precip_frac_2, wm,
            hm_mean, w_tol, hm_tol, a_exp, b_exp))

        # Calculate the covariance (overall) of two hydrometeors, <hmx'hmy'>,
        # for each unique set of two different hydrometeors.
        for hmy_idx in range(hm_idx + 1, hydromet_dim):

            # Unpack the mean (in-precip) of the second hydrometeor, hmy, in
            # each PDF component.
            mu_hmy_1 = hydromet_pdf_params.mu_hm_1[..., hmy_idx]
            mu_hmy_2 = hydromet_pdf_params.mu_hm_2[..., hmy_idx]

            # Unpack the standard deviation (in-precip) of hmy in each PDF
            # component.
            sigma_hmy_1 = hydromet_pdf_params.sigma_hm_1[..., hmy_idx]
            sigma_hmy_2 = hydromet_pdf_params.sigma_hm_2[..., hmy_idx]

            # Unpack the correlation (in-precip) of hm and hmy in each PDF
            # component.
            corr_hm_hmy_1 = hydromet_pdf_params.corr_hmx_hmy_1[..., hm_idx,hmy_idx]
            corr_hm_hmy_2 = hydromet_pdf_params.corr_hmx_hmy_2[..., hm_idx,hmy_idx]

            # Unpack the mean (overall) value of hmy.
            hmy_mean = hydromet[..., hmy_idx]

            # Unpack the tolerance value for the second hydrometeor, hmy.
            hmy_tol = hm_metadata.hydromet_tol[hmy_idx]

            # Calculate the covariance <hmx'hmy'>.
            hmxphmyp_zt = hmxphmyp_zt.at[..., hmy_idx,hm_idx].set(hmxphmyp_integral_covar(
                mu_hm_1, mu_hm_2, mu_hmy_1, mu_hmy_2,
                sigma_hm_1, sigma_hm_2, sigma_hmy_1, sigma_hmy_2,
                corr_hm_hmy_1, corr_hm_hmy_2, mixt_frac, precip_frac_1,
                precip_frac_2, hm_mean, hmy_mean, hm_tol,
                hmy_tol))


    # Statistics
    for hm_idx in range(hydromet_dim):
        hm_type = hm_metadata.hydromet_list[hm_idx]
        stats = stats.update("wp2" + hm_type[:2] + "p", wp2hmp[..., hm_idx])
        stats = stats.update("rtp" + hm_type[:2] + "p", zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, rtphmp_zt[..., hm_idx]))
        stats = stats.update("thlp" + hm_type[:2] + "p", zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, thlphmp_zt[..., hm_idx]))
    for hm_idx in range(hydromet_dim):
        for hmy_idx in range(hm_idx + 1, hydromet_dim):
            var_name = hm_metadata.hydromet_list[hm_idx][:2] + "p" + hm_metadata.hydromet_list[hmy_idx][:2] + "p"
            stats = stats.update(var_name, zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hmxphmyp_zt[..., hmy_idx, hm_idx]))
    return rtphmp_zt, thlphmp_zt, wp2hmp, stats


def xphmp_integral_covar(mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, sigma_x_1, sigma_x_2, sigma_hm_1, sigma_hm_2,
                         corr_x_hm_1, corr_x_hm_2, mixt_frac, precip_frac_1, precip_frac_2,
                         x_mean, x_tol, hm_tol):
    """Covariance <x'hm'> (mixed_moment_PDF_integrals.F90:xphmp_integral_covar) — the streamlined a=b=1 case.

    x is binormal (w/rt/thl/sclr), hm a precipitating hydrometeor. Within component i (in-precip only, since
    out-of-precip hm=0): E[(x-<x>)hm] = (μ_x_i - <x>) μ_hm_i + corr σ_x_i σ_hm_i. The Fortran's 4-way branch
    (drop the correlation term for whichever component has x or hm constant) decomposes into an independent
    per-component `jnp.where` — differentiable and equivalent (verified in the test). Uses hm's normal-space
    (in-precip) moments directly.
    """
    # Description:
    # Solves for the covariance < x'hm' >.  The variable "x" is a variable that
    # is distributed binormally (meaning that it has a normally-distributed
    # individual marginal in each PDF component).  This applies to such
    # variables as w, rt, thl, sclr, etc.  The variable "hm" stands for any
    # precipitating hydrometeor.  This applies to such variables as rr, Nr, ri,
    # Ni, etc.
    #
    # The covariance < x'hm' > can also be found by passing a = 1 and b = 1 into
    # function xp_a_hmpb_integrals_all_MM, which is found below.  However,
    # the code found here has been streamlined for the special case of
    # covariances.
    # References:
    #-----------------------------------------------------------------------
    varies_1 = (sigma_x_1 > x_tol) & (sigma_hm_1 > hm_tol)
    varies_2 = (sigma_x_2 > x_tol) & (sigma_hm_2 > hm_tol)
    return (mixt_frac * precip_frac_1 * ((mu_x_1 - x_mean) * mu_hm_1
                + jnp.where(varies_1, corr_x_hm_1 * sigma_x_1 * sigma_hm_1, 0.0))
            + (1.0 - mixt_frac) * precip_frac_2 * ((mu_x_2 - x_mean) * mu_hm_2
                + jnp.where(varies_2, corr_x_hm_2 * sigma_x_2 * sigma_hm_2, 0.0)))


def hmxphmyp_integral_covar(mu_hmx_1, mu_hmx_2, mu_hmy_1, mu_hmy_2, sigma_hmx_1, sigma_hmx_2,
                            sigma_hmy_1, sigma_hmy_2, corr_hmx_hmy_1, corr_hmx_hmy_2, mixt_frac,
                            precip_frac_1, precip_frac_2, hmx_mean, hmy_mean, hmx_tol, hmy_tol):
    """Covariance <hmx'hmy'> of two precipitating hydrometeors
    (mixed_moment_PDF_integrals.F90:hmxphmyp_integral_covar).

    E[hmx·hmy] over the mixture (in-precip only) minus <hmx><hmy>: within component i in-precip,
    E[hmx hmy] = μ_hmx_i μ_hmy_i + corr σ_hmx_i σ_hmy_i. Same per-component `jnp.where` decomposition of the
    Fortran 4-way branch.
    """
    # Description:
    # Solves for the covariance of two precipitating hydrometoers < hmx'hmy' >.
    # The variable "hmx" stands for any precipitating hydrometeor.  This applies
    # to such variables as rr, Nr, ri, Ni, etc.  The variable "hmy" stands for
    # any different precipitating hydrometeor.
    # References:
    #-----------------------------------------------------------------------
    varies_1 = (sigma_hmx_1 > hmx_tol) & (sigma_hmy_1 > hmy_tol)
    varies_2 = (sigma_hmx_2 > hmx_tol) & (sigma_hmy_2 > hmy_tol)
    return (mixt_frac * precip_frac_1 * (mu_hmx_1 * mu_hmy_1
                + jnp.where(varies_1, corr_hmx_hmy_1 * sigma_hmx_1 * sigma_hmy_1, 0.0))
            + (1.0 - mixt_frac) * precip_frac_2 * (mu_hmx_2 * mu_hmy_2
                + jnp.where(varies_2, corr_hmx_hmy_2 * sigma_hmx_2 * sigma_hmy_2, 0.0))
            - hmx_mean * hmy_mean)


def xp_a_hmpb_integrals_all_MM(mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, mu_hm_1_n, mu_hm_2_n,
                               sigma_x_1, sigma_x_2, sigma_hm_1, sigma_hm_2, sigma_hm_1_n, sigma_hm_2_n,
                               corr_x_hm_1_n, corr_x_hm_2_n, mixt_frac, precip_frac_1, precip_frac_2,
                               x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp):
    """Any mixed moment <x'^a hm'^b> over the full 2-component PDF
    (mixed_moment_PDF_integrals.F90:xp_a_hmpb_integrals_all_MM). x binormal, hm a precipitating hydrometeor."""
    # Description:
    # Solves for any mixed moment < x'^a hm'^b >, where a and b are both
    # integers with values greater than or equal to 0.  The variable "x" is a
    # variable that is distributed binormally (meaning that it has a
    # normally-distributed individual marginal in each PDF component).  This
    # applies to such variables as w, rt, thl, sclr, etc.  The variable "hm"
    # stands for any precipitating hydrometeor.  This applies to such variables
    # as rr, Nr, ri, Ni, etc.
    # References:
    #-----------------------------------------------------------------------
    return mixt_frac * bivar_NL_x_hm_all_MM_comp_eq(
        mu_x_1, mu_hm_1, mu_hm_1_n, sigma_x_1, sigma_hm_1, sigma_hm_1_n,
        corr_x_hm_1_n, precip_frac_1, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp) \
        + (1.0 - mixt_frac) * bivar_NL_x_hm_all_MM_comp_eq(
            mu_x_2, mu_hm_2, mu_hm_2_n, sigma_x_2, sigma_hm_2, sigma_hm_2_n,
            corr_x_hm_2_n, precip_frac_2, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp)


def bivar_NL_x_hm_all_MM_comp_eq(mu_x_i, mu_hm_i, mu_hm_i_n, sigma_x_i, sigma_hm_i, sigma_hm_i_n,
                                 corr_x_hm_i_n, precip_frac_i, x_mean, hm_mean, x_tol, hm_tol, a_exp, b_exp):
    """Per-PDF-component contribution to <x'^a hm'^b> (mixed_moment_PDF_integrals.F90:bivar_NL_x_hm_all_MM_comp_eq).

    x is normal-marginal (w/rt/thl/sclr...), hm is a precipitating hydrometeor (lognormal-marginal in-precip).
    The component splits into within-precip (fraction precip_frac_i, where hm follows its lognormal) and
    out-of-precip (1-precip_frac_i, where hm=0 → (hm-<hm>)^b = (-<hm>)^b, x still normal). The Fortran selects
    one of four closed forms depending on whether x and/or hm vary (σ vs tol); here all four are evaluated and
    selected with `jnp.where` (each is finite for finite inputs → safe value+gradient), so the port is
    vectorizable and differentiable while reproducing the exact branch the Fortran takes.

    Args mirror the Fortran; a_exp/b_exp are static non-negative integers.
    """
    # Description:
    # Solves the portion of the integral for < x'^a hm'^b > that relates to the
    # ith PDF component.  This takes into account the portions of the component
    # that are within-precipitation and outside-precipitation.  The equation
    # and function that are used depends on whether or not x and/or hm vary
    # inside of a PDF component.
    # References:
    #-----------------------------------------------------------------------
    one = 1.0
    out_precip = (one - precip_frac_i) * (-hm_mean) ** b_exp     # hm = 0 contribution (common to most branches)

    # (1) both x and hm constant in the component
    both_const = (mu_x_i - x_mean) ** a_exp * (
        precip_frac_i * (mu_hm_i - hm_mean) ** b_exp + out_precip)
    # (2) only x constant
    x_const = (mu_x_i - x_mean) ** a_exp * (
        precip_frac_i * univar_L_int_PDF_comp_all_MM(mu_hm_i_n, sigma_hm_i_n, hm_mean, b_exp) + out_precip)
    # (3) only hm constant
    univar_N_x = univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp)
    hm_const = (precip_frac_i * (mu_hm_i - hm_mean) ** b_exp + out_precip) * univar_N_x
    # (4) both vary
    both_vary = precip_frac_i * bivar_NL_int_PDF_comp_all_MM(
        mu_x_i, mu_hm_i_n, sigma_x_i, sigma_hm_i_n, corr_x_hm_i_n, x_mean, hm_mean, a_exp, b_exp) \
        + out_precip * univar_N_x

    x_is_const = sigma_x_i <= x_tol
    hm_is_const = sigma_hm_i <= hm_tol
    return jnp.where(x_is_const & hm_is_const, both_const,
                     jnp.where(x_is_const, x_const,
                               jnp.where(hm_is_const, hm_const, both_vary)))


def bivar_NL_int_PDF_comp_all_MM(mu_x1_i, mu_x2_i_n, sigma_x1_i, sigma_x2_i_n, corr_x1_x2_i_n,
                                 x1_mean, x2_mean, a_exp, b_exp):
    """Bivariate normal-lognormal central mixed moment within a PDF component.

    Evaluates  INT INT (x1-<x1>)^a (x2-<x2>)^b P_NL_i(x1,x2) dx2 dx1, where x1 has a normal marginal and x2 a
    lognormal marginal (in the component), correlated through corr(x1, ln x2):
      = SUM(p=0:floor(a/2)) SUM(q=0:b) [a!/((a-2p)!p!)] [b!/((b-q)!q!)]
          (½ σ_x1²)^p (μ_x1 - <x1> + ρ σ_x1 σ_x2_n q)^(a-2p)
          (-<x2>)^(b-q) exp( μ_x2_n q + ½ σ_x2_n² q² )
    (mixed_moment_PDF_integrals.F90:bivar_NL_int_PDF_comp_all_MM). Equivalent to exponential-tilting the joint
    normal by exp(q·ln x2): the lognormal raw moment factors out and x1's mean shifts by ρ σ_x1 σ_x2_n q.

    Args:
        mu_x1_i:        Mean of x1 in the component              [x1 units].
        mu_x2_i_n:      Mean of ln x2 in the component           [ln(x2 units)].
        sigma_x1_i:     Std dev of x1 in the component           [x1 units].
        sigma_x2_i_n:   Std dev of ln x2 in the component        [-].
        corr_x1_x2_i_n: Correlation of x1 and ln x2              [-].
        x1_mean:        Overall mean <x1>                        [x1 units].
        x2_mean:        Overall mean <x2>                        [x2 units].
        a_exp, b_exp:   Non-negative integer moment orders (static).
    """
    # Description:
    # This function is the evaluated form of the following integral:
    #
    # INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
    #                          * P_NL_i( x1, x2 ) dx2 dx1;
    #
    # where P_NL_i( x1, x2 ) is the functional form of a bivariate
    # normal-lognormal PDF (in the PDF component), x1 is a variable that has an
    # individual marginal that is distributed normally (in the PDF component),
    # and x2 is a variable that has an individual marginal that is distributed
    # lognormally (in the PDF component).  Additionally, <x1> is the overall
    # mean of x1, <x2> is the overall mean of x2, "a" is the integer (>= 0)
    # order of the mixed moment with respect to x1, and b is the integer (>= 0)
    # order of the mixed moment with respect to x2.
    #
    # When the integral is evaluated, the equation is:
    #
    # INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
    #                          * P_NL_i( x1, x2 ) dx2 dx1
    # = SUM( p = 0:floor(a/2) ) SUM ( q = 0:b )
    #   ( a! / ( ( a - 2p )! p! ) ) * ( b! / ( ( b - q )! q! ) )
    #   * [ (1/2) * sigma_x1_i**2 ]**p
    #   * ( mu_x1_i - <x1>
    #       + corr_x1_x2_i_n * sigma_x1_i * sigma_x2_i_n * q )**(a-2p)
    #   * ( -<x2> )**(b-q)
    #   * exp{ mu_x2_i_n * q + (1/2) * sigma_x2_i_n**2 * q**2 }.
    # References:
    #-----------------------------------------------------------------------
    total = 0.0
    for p in range(a_exp // 2 + 1):
        fac_p = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        for q in range(b_exp + 1):
            fac_q = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
            shifted = mu_x1_i - x1_mean + corr_x1_x2_i_n * sigma_x1_i * sigma_x2_i_n * q
            total = total + (fac_p * fac_q) * (0.5 * sigma_x1_i ** 2) ** p \
                * shifted ** (a_exp - 2 * p) \
                * (-x2_mean) ** (b_exp - q) \
                * jnp.exp(mu_x2_i_n * q + 0.5 * sigma_x2_i_n ** 2 * q ** 2)
    return total


def univar_N_int_PDF_comp_all_MM(mu_x_i, sigma_x_i, x_mean, a_exp):
    """Central moment of a normally-distributed variable within a PDF component.

    Evaluates  INT(-inf:inf) (x - <x>)^a P_N_i(x) dx
      = SUM(p=0:floor(a/2)) [ a! / ((a-2p)! p!) ] (½ σ_i²)^p (μ_i - <x>)^(a-2p)
    (mixed_moment_PDF_integrals.F90:univar_N_int_PDF_comp_all_MM).

    Args:
        mu_x_i:    Mean of x in the component       [x units].
        sigma_x_i: Std dev of x in the component    [x units].
        x_mean:    Overall mean <x>                 [x units].
        a_exp:     Non-negative integer moment order (static).
    """
    # Description:
    # This function is the evaluated form of the following integral:
    #
    # INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx;
    #
    # where P_N_i( x ) is the functional form of a single-variable normal PDF
    # (in the PDF component), and x is a variable that has an individual
    # marginal that is distributed normally (in the PDF component).
    # Additionally, <x> is the overall mean of x, and "a" is the integer (>= 0)
    # order of the mixed moment with respect to x.
    #
    # When the integral is evaluated, the equation is:
    #
    # INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx
    # = SUM( p = 0:floor(a/2) ) ( a! / ( ( a - 2p )! p! ) )
    #   * [ (1/2) * sigma_x_i**2 ]**p * ( mu_x_i - <x> )**(a-2p).
    # References:
    #-----------------------------------------------------------------------
    total = 0.0
    for p in range(a_exp // 2 + 1):
        fac = math.factorial(a_exp) // (math.factorial(a_exp - 2 * p) * math.factorial(p))
        total = total + fac * (0.5 * sigma_x_i ** 2) ** p * (mu_x_i - x_mean) ** (a_exp - 2 * p)
    return total


def univar_L_int_PDF_comp_all_MM(mu_x_i_n, sigma_x_i_n, x_mean, b_exp):
    """Central moment of a lognormally-distributed variable within a PDF component.

    Evaluates  INT(0:inf) (x - <x>)^b P_L_i(x) dx
      = SUM(q=0:b) [ b! / ((b-q)! q!) ] (-<x>)^(b-q) exp( μ_n q + ½ σ_n² q² )
    where (μ_n, σ_n) are the mean/std of ln x in the component
    (mixed_moment_PDF_integrals.F90:univar_L_int_PDF_comp_all_MM).

    Args:
        mu_x_i_n:    Mean of ln x in the component   [ln(x units)].
        sigma_x_i_n: Std dev of ln x in the component [-].
        x_mean:      Overall mean <x>                [x units].
        b_exp:       Non-negative integer moment order (static).
    """
    # Description:
    # This function is the evaluated form of the following integral:
    #
    # INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
    #
    # where P_L_i( x ) is the functional form of a single-variable lognormal PDF
    # (in the PDF component), and x is a variable that has an individual
    # marginal that is distributed lognormally (in the PDF component).
    # Additionally, <x> is the overall mean of x, and b is the integer (>= 0)
    # order of the mixed moment with respect to x.
    #
    # When the integral is evaluated, the equation is:
    #
    # INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
    # = SUM ( q = 0:b ) ( b! / ( ( b - q )! q! ) )
    #   * ( -<x> )**(b-q)
    #   * exp{ mu_x_i_n * q + (1/2) * sigma_x_i_n**2 * q**2 }.
    # References:
    #-----------------------------------------------------------------------
    total = 0.0
    for q in range(b_exp + 1):
        fac = math.factorial(b_exp) // (math.factorial(b_exp - q) * math.factorial(q))
        total = total + fac * (-x_mean) ** (b_exp - q) * jnp.exp(
            mu_x_i_n * q + 0.5 * sigma_x_i_n ** 2 * q ** 2)
    return total
