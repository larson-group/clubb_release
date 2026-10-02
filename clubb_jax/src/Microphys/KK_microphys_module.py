"""KK microphysics interface, mirroring KK_microphys_module.F90.

Columns and levels are batched in JAX; output/inout arguments are returned.
The numerical kernels in KK_microphys remain independently implemented.
"""
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Cp, Rd, Rv, rho_lw, pi, rr_tol, Nr_tol, eps, cm3_per_m3
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq
from clubb_jax.src.Microphys import parameters_microphys
from clubb_jax.src.Microphys.KK_microphys import parameters_KK
from clubb_jax.src.Microphys.KK_microphys.KK_utilities import G_T_p

def KK_local_microphys_driver(gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata,
                              l_latin_hypercube, thlm, wm_zt, p_in_Pa, exner, rho,
                              cloud_frac, w_std_dev, dzq, rcm, Ncm, chi, rvm,
                              hydromet, saturation_formula, sample_weight, stats):
    # Run local K-K microphysics for all model columns and update statistics
    # with full-column fields after the microphysics calculations finish.
    # SILHS remains rejected at initialization; its dummy driver must be ported
    # before weighted/noninteractive sample statistics can be enabled.
    if l_latin_hypercube:
        raise NotImplementedError('SILHS microphysics is disabled')
    (hydromet_mc, hydromet_vel, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
     rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag,
     mvrr_diag, rrm_src_adj_diag, Nrm_src_adj_diag, rrm_evap_adj_diag,
     Nrm_evap_adj_diag, rrm_mc_nonadj_diag) = KK_local_microphys_core(
        gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata, l_latin_hypercube,
        thlm, p_in_Pa, exner, rho, rcm, Ncm, chi, hydromet, saturation_formula)

    # Fortran internal procedure: Python must define it before the first call.
    # The closure carries immutable JaxStats in place of host association.
    def update_microphys_stat(name, value):
        nonlocal stats
        stats = stats.update(name, value)

    # Output values for statistics.
    update_microphys_stat('rrm_evap', rrm_evap_diag)
    update_microphys_stat('rrm_auto', rrm_auto_diag)
    update_microphys_stat('rrm_accr', rrm_accr_diag)
    update_microphys_stat('mvrr', mvrr_diag)
    update_microphys_stat('Nrm_evap', Nrm_evap_diag)
    update_microphys_stat('Nrm_auto', Nrm_auto_diag)
    update_microphys_stat('rrm_src_adj', rrm_src_adj_diag)
    update_microphys_stat('Nrm_src_adj', Nrm_src_adj_diag)
    update_microphys_stat('rrm_evap_adj', rrm_evap_adj_diag)
    update_microphys_stat('Nrm_evap_adj', Nrm_evap_adj_diag)
    update_microphys_stat('rrm_mc_nonadj', rrm_mc_nonadj_diag)
    return (stats, hydromet_mc, hydromet_vel, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
            rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag)


def KK_local_microphys_core(gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata,
                            l_latin_hypercube, thlm, p_in_Pa, exner, rho,
                            rcm, Ncm, chi, hydromet, saturation_formula):
    # Description:
    # References:
    # Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    #    Parameterization in a Large-Eddy Simulation Model of Marine
    #    Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    #-----------------------------------------------------------------------
    from clubb_jax.src.Microphys.KK_microphys.KK_local_means import (
        KK_mvr_local_mean, KK_evap_local_mean, KK_auto_local_mean, KK_accr_local_mean)
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import KK_Nrm_evap_local_mean, KK_Nrm_auto_mean
    from clubb_jax.src.Microphys.advance_microphys_module import get_cloud_top_level
    from clubb_jax.src.CLUBB_core.constants_clubb import Nc_tol
    # Set up mean fields and microphysics tendency adjustment flags.
    rrm = hydromet[..., hm_metadata.iirr]
    Nrm = hydromet[..., hm_metadata.iiNr]
    l_src_adj_enabled = l_evap_adj_enabled = l_clip_positive_sed = True
    Ncm_mc = jnp.zeros_like(rcm)
    if parameters_microphys.l_silhs_KK_convergence_adj_mean and l_latin_hypercube:
        l_src_adj_enabled = l_evap_adj_enabled = l_clip_positive_sed = False
    # Microphysics tendency loop, batched over levels and columns.
    KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef = KK_tendency_coefs(thlm, exner, p_in_Pa, rho, saturation_formula)
    KK_mean_vol_rad = jnp.where(rrm > rr_tol, KK_mvr_local_mean(rrm, Nrm, KK_mvr_coef, Nr_tol), 0.0)
    KK_evap_tndcy = jnp.where((rrm > rr_tol) & (Nrm > Nr_tol), KK_evap_local_mean(chi, rrm, Nrm, KK_evap_coef), 0.0)
    KK_auto_tndcy = jnp.where(Ncm > Nc_tol, KK_auto_local_mean(chi, jnp.maximum(Ncm, Nc_tol), KK_auto_coef), 0.0)
    KK_accr_tndcy = jnp.where(rrm > rr_tol, KK_accr_local_mean(chi, rrm, KK_accr_coef), 0.0)
    KK_Nrm_evap_tndcy = jnp.where((rrm > rr_tol) & (Nrm > Nr_tol), KK_Nrm_evap_local_mean(KK_evap_tndcy, Nrm, rrm, dt), 0.0)
    KK_Nrm_auto_tndcy = KK_Nrm_auto_mean(KK_auto_tndcy)
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, adj_terms = KK_microphys_adjust(
        dt, exner, rcm, rrm, Nrm, KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy,
        KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy, l_src_adj_enabled, l_evap_adj_enabled)
    # Boundary conditions for microphysics tendencies.
    rrm_mc, Nrm_mc = rrm_mc.at[:, -1].set(0.0), Nrm_mc.at[:, -1].set(0.0)
    KK_mean_vol_rad = KK_mean_vol_rad.at[:, -1].set(0.0)
    rvm_mc, rcm_mc, thlm_mc = rvm_mc.at[:, -1].set(0.0), rcm_mc.at[:, -1].set(0.0), thlm_mc.at[:, -1].set(0.0)
    cloud_top_level = get_cloud_top_level(nzt, ngrdcol, rcm, hydromet, hydromet_dim,
            hm_metadata.iiri)
    Vrr, VNr = KK_upscaled_sedimentation(ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad, l_clip_positive_sed)
    hydromet_mc, hydromet_vel = KK_microphys_output(ngrdcol, nzt, hydromet_dim, hm_metadata, Vrr, VNr, rrm_mc, Nrm_mc)
    rrm_auto_diag, rrm_accr_diag, rrm_evap_diag = KK_auto_tndcy, KK_accr_tndcy, KK_evap_tndcy
    Nrm_auto_diag, Nrm_evap_diag = KK_Nrm_auto_tndcy, KK_Nrm_evap_tndcy
    mvrr_diag = KK_mean_vol_rad
    rrm_src_adj_diag, Nrm_src_adj_diag, rrm_evap_adj_diag, Nrm_evap_adj_diag = adj_terms
    rrm_mc_nonadj_diag = KK_auto_tndcy + KK_accr_tndcy + KK_evap_tndcy
    return (hydromet_mc, hydromet_vel, Ncm_mc, rcm_mc, rvm_mc, thlm_mc,
            rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, Nrm_auto_diag, Nrm_evap_diag,
            mvrr_diag, rrm_src_adj_diag, Nrm_src_adj_diag, rrm_evap_adj_diag,
            Nrm_evap_adj_diag, rrm_mc_nonadj_diag)


def KK_upscaled_microphys(gr, ngrdcol, dt, nzt, nzm, pdf_dim, hydromet_dim, hm_metadata,
                          wm_zt, rtm, thlm, p_in_Pa, exner, rho, rcm,
                          pdf_params, hydromet_pdf_params, precip_fracs, hydromet,
                          mu_x_1_n, mu_x_2_n, sigma_x_1_n, sigma_x_2_n,
                          corr_array_1_n, corr_array_2_n, saturation_formula, stats):
    # Description:
    # Version of KK microphysics scheme that is analytically upscaled by
    # integrating over the product of the microphysics tendency and the
    # functional form of the PDF.
    # References:
    # Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of a local
    #    microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol. Soc.,
    #    139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
    #
    # Griffin, B. M. and V. E. Larson, 2013:  Analytic upscaling of a local
    #    microphysics scheme. Part II: Simulations.  Q. J. Roy. Meteorol. Soc.,
    #    139, 670, 58--69, doi:http://dx.doi.org/10.1002/qj.1966.
    #
    # Griffin, B. M., 2016:  Improving the Subgrid-Scale Representation of
    #    Hydrometeors and Microphysical Feedback Effects Using a Multivariate
    #    PDF.  Doctoral dissertation, University of Wisconsin -- Milwaukee,
    #    Milwaukee, WI, Paper 1144, 165 pp., URL
    #    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    #
    # Griffin, B. M. and V. E. Larson, 2016:  Supplement of A new subgrid-scale
    #    representation of hydrometeor fields using a multivariate PDF.
    #    Geosci. Model Dev., 9, 6,
    #    doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    #
    # Griffin, B. M. and V. E. Larson, 2016:  Parameterizing microphysical
    #    effects on variances and covariances of moisture and heat content using
    #    a multivariate probability density function: a study with CLUBB (tag
    #    MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    #    doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    #-----------------------------------------------------------------------
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_means import (
        KK_evap_upscaled_mean, KK_auto_upscaled_mean, KK_accr_upscaled_mean, KK_mvr_upscaled_mean)
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import KK_upscaled_covar_driver
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import KK_sed_vel_covars
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import KK_Nrm_evap_upscaled_mean, KK_Nrm_auto_mean
    from clubb_jax.src.Microphys.advance_microphys_module import get_cloud_top_level

    rrm = hydromet[..., hm_metadata.iirr]
    Nrm = hydromet[..., hm_metadata.iiNr]
    l_src_adj_enabled = l_evap_adj_enabled = True
    # Setup mixture fraction.
    mixt_frac = pdf_params.mixt_frac
    # Microphysics tendency loop: batched over columns and thermodynamic levels.
    KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef = KK_tendency_coefs(
        thlm, exner, p_in_Pa, rho, saturation_formula)
    # Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
    mu_w_1        = mu_x_1_n[..., hm_metadata.iiPDF_w]
    mu_w_2        = mu_x_2_n[..., hm_metadata.iiPDF_w]
    mu_chi_1      = mu_x_1_n[..., hm_metadata.iiPDF_chi]
    mu_chi_2      = mu_x_2_n[..., hm_metadata.iiPDF_chi]
    mu_eta_1      = mu_x_1_n[..., hm_metadata.iiPDF_eta]
    mu_eta_2      = mu_x_2_n[..., hm_metadata.iiPDF_eta]
    mu_rr_1_n     = mu_x_1_n[..., hm_metadata.iiPDF_rr]
    mu_rr_2_n     = mu_x_2_n[..., hm_metadata.iiPDF_rr]
    mu_Nr_1_n     = mu_x_1_n[..., hm_metadata.iiPDF_Nr]
    mu_Nr_2_n     = mu_x_2_n[..., hm_metadata.iiPDF_Nr]
    mu_Ncn_1_n    = mu_x_1_n[..., hm_metadata.iiPDF_Ncn]
    mu_Ncn_2_n    = mu_x_2_n[..., hm_metadata.iiPDF_Ncn]
    sigma_w_1     = sigma_x_1_n[..., hm_metadata.iiPDF_w]
    sigma_w_2     = sigma_x_2_n[..., hm_metadata.iiPDF_w]
    sigma_chi_1   = sigma_x_1_n[..., hm_metadata.iiPDF_chi]
    sigma_chi_2   = sigma_x_2_n[..., hm_metadata.iiPDF_chi]
    sigma_eta_1   = sigma_x_1_n[..., hm_metadata.iiPDF_eta]
    sigma_eta_2   = sigma_x_2_n[..., hm_metadata.iiPDF_eta]
    sigma_rr_1_n  = sigma_x_1_n[..., hm_metadata.iiPDF_rr]
    sigma_rr_2_n  = sigma_x_2_n[..., hm_metadata.iiPDF_rr]
    sigma_Nr_1_n  = sigma_x_1_n[..., hm_metadata.iiPDF_Nr]
    sigma_Nr_2_n  = sigma_x_2_n[..., hm_metadata.iiPDF_Nr]
    sigma_Ncn_1_n = sigma_x_1_n[..., hm_metadata.iiPDF_Ncn]
    sigma_Ncn_2_n = sigma_x_2_n[..., hm_metadata.iiPDF_Ncn]

    # Unpack variables from hydromet_pdf_params
    mu_rr_1     = hydromet_pdf_params.mu_hm_1[..., hm_metadata.iirr]
    mu_rr_2     = hydromet_pdf_params.mu_hm_2[..., hm_metadata.iirr]
    mu_Nr_1     = hydromet_pdf_params.mu_hm_1[..., hm_metadata.iiNr]
    mu_Nr_2     = hydromet_pdf_params.mu_hm_2[..., hm_metadata.iiNr]
    mu_Ncn_1    = hydromet_pdf_params.mu_Ncn_1
    mu_Ncn_2    = hydromet_pdf_params.mu_Ncn_2
    sigma_rr_1  = hydromet_pdf_params.sigma_hm_1[..., hm_metadata.iirr]
    sigma_rr_2  = hydromet_pdf_params.sigma_hm_2[..., hm_metadata.iirr]
    sigma_Nr_1  = hydromet_pdf_params.sigma_hm_1[..., hm_metadata.iiNr]
    sigma_Nr_2  = hydromet_pdf_params.sigma_hm_2[..., hm_metadata.iiNr]
    sigma_Ncn_1 = hydromet_pdf_params.sigma_Ncn_1
    sigma_Ncn_2 = hydromet_pdf_params.sigma_Ncn_2

    rr_1          = hydromet_pdf_params.hm_1[..., hm_metadata.iirr]
    rr_2          = hydromet_pdf_params.hm_2[..., hm_metadata.iirr]
    Nr_1          = hydromet_pdf_params.hm_1[..., hm_metadata.iiNr]
    Nr_2          = hydromet_pdf_params.hm_2[..., hm_metadata.iiNr]

    # Unpack corr_array_1_n into correlations (1st PDF component).
    corr_chi_eta_1   = corr_array_1_n[..., hm_metadata.iiPDF_eta, hm_metadata.iiPDF_chi]
    corr_w_chi_1     = corr_array_1_n[..., hm_metadata.iiPDF_w,hm_metadata.iiPDF_chi]
    corr_chi_rr_1_n  = corr_array_1_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_chi]
    corr_chi_Nr_1_n  = corr_array_1_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_chi]
    corr_chi_Ncn_1_n = corr_array_1_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_chi]
    corr_eta_rr_1_n  = corr_array_1_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_eta]
    corr_eta_Nr_1_n  = corr_array_1_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_eta]
    corr_eta_Ncn_1_n = corr_array_1_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_eta]
    corr_w_rr_1_n    = corr_array_1_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_w]
    corr_w_Nr_1_n    = corr_array_1_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_w]
    corr_w_Ncn_1_n   = corr_array_1_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_w]
    corr_rr_Nr_1_n   = corr_array_1_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_rr]

    # Unpack corr_array_2_n into correlations (2nd PDF component).
    corr_chi_eta_2   = corr_array_2_n[..., hm_metadata.iiPDF_eta, hm_metadata.iiPDF_chi]
    corr_w_chi_2     = corr_array_2_n[..., hm_metadata.iiPDF_w,hm_metadata.iiPDF_chi]
    corr_chi_rr_2_n  = corr_array_2_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_chi]
    corr_chi_Nr_2_n  = corr_array_2_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_chi]
    corr_chi_Ncn_2_n = corr_array_2_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_chi]
    corr_eta_rr_2_n  = corr_array_2_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_eta]
    corr_eta_Nr_2_n  = corr_array_2_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_eta]
    corr_eta_Ncn_2_n = corr_array_2_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_eta]
    corr_w_rr_2_n    = corr_array_2_n[..., hm_metadata.iiPDF_rr, hm_metadata.iiPDF_w]
    corr_w_Nr_2_n    = corr_array_2_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_w]
    corr_w_Ncn_2_n   = corr_array_2_n[..., hm_metadata.iiPDF_Ncn, hm_metadata.iiPDF_w]
    corr_rr_Nr_2_n   = corr_array_2_n[..., hm_metadata.iiPDF_Nr, hm_metadata.iiPDF_rr]


    precip_frac_1 = precip_fracs.precip_frac_1
    precip_frac_2 = precip_fracs.precip_frac_2
    # Calculate the values of the upscaled KK microphysics tendencies.
    KK_evap_tndcy = KK_evap_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1,
        mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1,
        sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n,
        corr_rr_Nr_2_n, KK_evap_coef, mixt_frac, precip_frac_1, precip_frac_2)
    KK_auto_tndcy = KK_auto_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n,
        mu_Ncn_2_n, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2,
        sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, KK_auto_coef,
        mixt_frac)
    KK_accr_tndcy = KK_accr_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n,
        mu_rr_2_n, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2,
        sigma_rr_1_n, sigma_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, mixt_frac,
        precip_frac_1, precip_frac_2)
    KK_mean_vol_rad = KK_mvr_upscaled_mean(
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
        mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, sigma_rr_2,
        sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac, precip_frac_1,
        precip_frac_2)
    # Calculate the implicit and explicit turbulent sedimentation terms.
    sedimentation = KK_sed_vel_covars(
        rrm, rr_1, rr_2, Nrm, Nr_1, Nr_2, KK_mean_vol_rad,
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
        mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, sigma_rr_2,
        sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac)
    if parameters_microphys.l_var_covar_src:
        (wprtp_mc_zt, wpthlp_mc_zt, rtp2_mc_zt, thlp2_mc_zt, rtpthlp_mc_zt,
         w_KK_evap_covar, rt_KK_evap_covar, thl_KK_evap_covar,
         w_KK_auto_covar, rt_KK_auto_covar, thl_KK_auto_covar,
         w_KK_accr_covar, rt_KK_accr_covar, thl_KK_accr_covar) = KK_upscaled_covar_driver(
        wm_zt, rtm, thlm, exner, mu_w_1,
        mu_w_2, mu_chi_1, mu_chi_2, mu_eta_1, mu_eta_2,
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_Ncn_1,
        mu_Ncn_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        mu_Ncn_1_n, mu_Ncn_2_n, sigma_w_1, sigma_w_2, sigma_chi_1,
        sigma_chi_2, sigma_eta_1, sigma_eta_2, sigma_rr_1, sigma_rr_2,
        sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, sigma_rr_1_n,
        sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n,
        corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, corr_w_rr_2_n, corr_w_Nr_1_n,
        corr_w_Nr_2_n, corr_w_Ncn_1_n, corr_w_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2,
        corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_chi_Ncn_1_n,
        corr_chi_Ncn_2_n, corr_eta_rr_1_n, corr_eta_rr_2_n, corr_eta_Nr_1_n, corr_eta_Nr_2_n,
        corr_eta_Ncn_1_n, corr_eta_Ncn_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac,
        precip_frac_1, precip_frac_2, KK_evap_coef, KK_auto_coef, KK_accr_coef,
        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy, pdf_params.rt_1, pdf_params.rt_2,
        pdf_params.thl_1, pdf_params.thl_2, pdf_params.crt_1, pdf_params.crt_2, pdf_params.cthl_1,
        pdf_params.cthl_2)
    # KK rain drop concentration microphysics tendencies.
    KK_Nrm_evap_tndcy = KK_Nrm_evap_upscaled_mean(
        mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1,
        mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1,
        sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n,
        corr_rr_Nr_2_n, KK_evap_coef, mixt_frac, precip_frac_1, precip_frac_2,
        dt)
    KK_Nrm_auto_tndcy = KK_Nrm_auto_mean(KK_auto_tndcy)
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, adj_terms = KK_microphys_adjust(
        dt, exner, rcm, rrm, Nrm, KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy,
        KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy, l_src_adj_enabled, l_evap_adj_enabled)
    # Calculate the variance of KK rain drop mean volume radius.
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_variances import variance_KK_mvr
    if stats.l_sample and stats.var_on_stats_list("KK_mvr_variance_zt"):
        KK_mvr_variance = variance_KK_mvr(
            mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
            mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1, sigma_rr_2,
            sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
            sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mean_vol_rad,
            KK_mvr_coef, mixt_frac, precip_frac_1, precip_frac_2)
    # Statistics use full-column fields after all tendencies are calculated.
    stats = stats.update('rr_KK_mvr_covar_zt', sedimentation['rr_KK_mvr_covar'])
    stats = stats.update('Nr_KK_mvr_covar_zt', sedimentation['Nr_KK_mvr_covar'])
    if parameters_microphys.l_var_covar_src:
        stats = stats.update('w_KK_evap_covar_zt', w_KK_evap_covar)
        stats = stats.update('rt_KK_evap_covar_zt', rt_KK_evap_covar)
        stats = stats.update('thl_KK_evap_covar_zt', thl_KK_evap_covar)
        stats = stats.update('w_KK_auto_covar_zt', w_KK_auto_covar)
        stats = stats.update('rt_KK_auto_covar_zt', rt_KK_auto_covar)
        stats = stats.update('thl_KK_auto_covar_zt', thl_KK_auto_covar)
        stats = stats.update('w_KK_accr_covar_zt', w_KK_accr_covar)
        stats = stats.update('rt_KK_accr_covar_zt', rt_KK_accr_covar)
        stats = stats.update('thl_KK_accr_covar_zt', thl_KK_accr_covar)
    stats = stats.update('rrm_evap', KK_evap_tndcy)
    stats = stats.update('rrm_auto', KK_auto_tndcy)
    stats = stats.update('rrm_accr', KK_accr_tndcy)
    stats = stats.update('mvrr', KK_mean_vol_rad)
    stats = stats.update('Nrm_evap', KK_Nrm_evap_tndcy)
    stats = stats.update('Nrm_auto', KK_Nrm_auto_tndcy)
    if stats.l_sample and stats.var_on_stats_list("KK_mvr_variance_zt"):
        stats = stats.update('KK_mvr_variance_zt', KK_mvr_variance)
    if parameters_microphys.l_var_covar_src:
        wprtp_mc = zt2zm(nzm, nzt, gr.ngrdcol, gr, wprtp_mc_zt).at[:, 0].set(0.0).at[:, -1].set(0.0)
        wpthlp_mc = zt2zm(nzm, nzt, gr.ngrdcol, gr, wpthlp_mc_zt).at[:, 0].set(0.0).at[:, -1].set(0.0)
        rtp2_mc = zt2zm(nzm, nzt, gr.ngrdcol, gr, rtp2_mc_zt).at[:, 0].set(0.0).at[:, -1].set(0.0)
        thlp2_mc = zt2zm(nzm, nzt, gr.ngrdcol, gr, thlp2_mc_zt).at[:, 0].set(0.0).at[:, -1].set(0.0)
        rtpthlp_mc = zt2zm(nzm, nzt, gr.ngrdcol, gr, rtpthlp_mc_zt).at[:, 0].set(0.0).at[:, -1].set(0.0)
    else:
        wprtp_mc = wpthlp_mc = rtp2_mc = thlp2_mc = rtpthlp_mc = jnp.zeros((gr.ngrdcol, nzm))
    # Boundary conditions for microphysics tendencies.
    rrm_mc = rrm_mc.at[:, -1].set(0.0)
    Nrm_mc = Nrm_mc.at[:, -1].set(0.0)
    KK_mean_vol_rad = KK_mean_vol_rad.at[:, -1].set(0.0)
    rvm_mc = rvm_mc.at[:, -1].set(0.0)
    rcm_mc = rcm_mc.at[:, -1].set(0.0)
    thlm_mc = thlm_mc.at[:, -1].set(0.0)
    cloud_top_level = get_cloud_top_level(nzt, ngrdcol, rcm, hydromet, hydromet_dim,
            hm_metadata.iiri)
    Vrr, VNr = KK_upscaled_sedimentation(ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad, True)
    hydromet_mc = jnp.zeros_like(hydromet)
    hydromet_vel = jnp.zeros_like(hydromet)
    hydromet_mc = hydromet_mc.at[..., hm_metadata.iirr].set(rrm_mc)
    hydromet_mc = hydromet_mc.at[..., hm_metadata.iiNr].set(Nrm_mc)
    hydromet_vel = hydromet_vel.at[..., hm_metadata.iirr].set(Vrr)
    hydromet_vel = hydromet_vel.at[..., hm_metadata.iiNr].set(VNr)
    # Turbulent sedimentation above cloud top and through the model top is zero.
    above = (jnp.arange(nzt)[None, :] > cloud_top_level[:, None]) & (cloud_top_level[:, None] > 0)
    hydromet_vel_covar_zt_impc = jnp.stack((sedimentation['Vrrprrp_impc'], sedimentation['VNrpNrp_impc']), axis=-1)
    hydromet_vel_covar_zt_expc = jnp.stack((sedimentation['Vrrprrp_expc'], sedimentation['VNrpNrp_expc']), axis=-1)
    hydromet_vel_covar_zt_impc = jnp.where(above[..., None], 0.0, hydromet_vel_covar_zt_impc).at[:, -1, :].set(0.0)
    hydromet_vel_covar_zt_expc = jnp.where(above[..., None], 0.0, hydromet_vel_covar_zt_expc).at[:, -1, :].set(0.0)
    stats = stats.update('Vrrprrp_expcalc', zt2zm(nzm, nzt, gr.ngrdcol, gr,
        hydromet_vel_covar_zt_impc[..., hm_metadata.iirr] * rrm + hydromet_vel_covar_zt_expc[..., hm_metadata.iirr]))
    stats = stats.update('VNrpNrp_expcalc', zt2zm(nzm, nzt, gr.ngrdcol, gr,
        hydromet_vel_covar_zt_impc[..., hm_metadata.iiNr] * Nrm + hydromet_vel_covar_zt_expc[..., hm_metadata.iiNr]))
    stats = stats.update('rrm_src_adj', adj_terms[0])
    stats = stats.update('Nrm_src_adj', adj_terms[1])
    stats = stats.update('rrm_evap_adj', adj_terms[2])
    stats = stats.update('Nrm_evap_adj', adj_terms[3])
    stats = stats.update('rrm_mc_nonadj', KK_auto_tndcy + KK_accr_tndcy + KK_evap_tndcy)
    return (stats, hydromet_mc, hydromet_vel, rcm_mc, rvm_mc, thlm_mc,
            hydromet_vel_covar_zt_impc, hydromet_vel_covar_zt_expc,
            wprtp_mc, wpthlp_mc, rtp2_mc, thlp2_mc, rtpthlp_mc)

def KK_tendency_coefs(thlm, exner, p_in_Pa, rho, saturation_formula):
    # Liquid water temperature.
    # Description:
    # References:
    # Eq. (3), Eq. (22), Eq. (29), and Eq. (33) of Khairoutdinov, M. and
    # Y. Kogan, 2000:  A New Cloud Physics Parameterization in a Large-Eddy
    # Simulation Model of Marine Stratocumulus.  Mon. Wea. Rev., 128, 229--243.
    #
    # Eq. (22), Eq. (28), Eq. (38), and Eq. (51) of Larson, V. E. and
    # B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
    # Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
    # doi:http://dx.doi.org/10.1002/qj.1967.
    #
    # Eq. (C21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
    # Representation of Hydrometeors and Microphysical Feedback Effects Using a
    # Multivariate PDF.  Doctoral dissertation, University of
    # Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
    # http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
    #
    # Eq. (S21) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
    # A new subgrid-scale representation of hydrometeor fields using a
    # multivariate PDF.  Geosci. Model Dev., 9, 6,
    # doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
    #
    # Eq. (A27) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
    # microphysical effects on variances and covariances of moisture and heat
    # content using a multivariate probability density function: a study with
    # CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
    # doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
    #-----------------------------------------------------------------------
    T_liq_in_K = thlm * exner
    r_sl = sat_mixrat_liq(p_in_Pa, T_liq_in_K, saturation_formula)
    Beta_Tl = (Rd / Rv) * (Lv / (Rd * T_liq_in_K)) * (Lv / (Cp * T_liq_in_K))
    KK_evap_coef = (3.0 * parameters_KK.C_evap * G_T_p(T_liq_in_K, p_in_Pa, saturation_formula)
                   * ((4.0 / 3.0) * pi * rho_lw) ** (2.0 / 3.0) * ((1.0 + Beta_Tl * r_sl) / r_sl))
    KK_auto_coef = 1350.0 * (rho / cm3_per_m3) ** parameters_KK.KK_auto_Nc_exp
    KK_accr_coef = 67.0
    KK_mvr_coef = ((4.0 / 3.0) * pi * rho_lw) ** (-1.0 / 3.0)
    return KK_evap_coef, KK_auto_coef, KK_accr_coef, KK_mvr_coef

def KK_microphys_adjust(dt, exner, rcm, rrm, Nrm,
                        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy,
                        KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy,
                        l_src_adj_enabled, l_evap_adj_enabled):
    """Assemble the KK microphysics state tendencies from the process rates.
    KK_microphys_module.F90:1196 (the upscaled path enables both adjustments).

    Source adjustment: limit auto+accr so they don't draw more cloud water than available
    (rate <= rcm/dt). Evaporation adjustment: limit so rain can't go negative (>= -rrm/dt,
    -Nrm/dt). Returns (rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc)."""
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import (
        KK_Nrm_auto_mean, KK_Nrm_evap_local_mean)
    from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Cp

    rrm_src_adj = Nrm_src_adj = jnp.zeros_like(rrm)
    rrm_source = KK_auto_tndcy + KK_accr_tndcy
    Nrm_source = KK_Nrm_auto_tndcy

    if l_src_adj_enabled:
        # Over a long step auto+accr may over-deplete rcm; cap the total source at rcm/dt.
        over = (rrm_source * dt) > rcm
        rrm_src_max = rcm / dt
        src_safe = jnp.where(rrm_source != 0.0, rrm_source, 1.0)
        rrm_auto_ratio = KK_auto_tndcy / src_safe
        rrm_src_adj = rrm_src_max - rrm_source
        Nrm_src_adj = KK_Nrm_auto_mean(rrm_auto_ratio * rrm_src_adj)
        rrm_source = jnp.where(over, rrm_src_max, rrm_source)
        Nrm_source = jnp.where(over, Nrm_source + Nrm_src_adj, Nrm_source)
        rrm_src_adj = jnp.where(over, rrm_src_adj, 0.0)
        Nrm_src_adj = jnp.where(over, Nrm_src_adj, 0.0)

    if l_evap_adj_enabled:
        rrm_evap_net = jnp.maximum(KK_evap_tndcy, -rrm / dt)
        # recompute Nrm evap from the net rrm evap when the rrm evap was limited
        limited = (jnp.abs(KK_evap_tndcy - rrm_evap_net)
                   > jnp.abs(KK_evap_tndcy + rrm_evap_net) * eps / 2.0) \
                  & (rrm > rr_tol) & (Nrm > Nr_tol)
        Nrm_evap_recomp = KK_Nrm_evap_local_mean(rrm_evap_net, Nrm, rrm, dt)
        Nrm_evap_net = jnp.where(limited, Nrm_evap_recomp, KK_Nrm_evap_tndcy)
        Nrm_evap_net = jnp.maximum(Nrm_evap_net, -Nrm / dt)
    else:
        rrm_evap_net = KK_evap_tndcy
        Nrm_evap_net = KK_Nrm_evap_tndcy

    rrm_mc = rrm_evap_net + rrm_source
    Nrm_mc = Nrm_evap_net + Nrm_source
    rvm_mc = -rrm_evap_net
    rcm_mc = -rrm_source
    thlm_mc = (Lv / (Cp * exner)) * rrm_mc
    return rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, (rrm_src_adj, Nrm_src_adj, rrm_evap_net - KK_evap_tndcy, Nrm_evap_net - KK_Nrm_evap_tndcy)

def KK_upscaled_sedimentation(ngrdcol, nzt, cloud_top_level, KK_mean_vol_rad, l_clip_positive_sed):
    # Description:
    # References:
    # Eq. (37) of Khairoutdinov, M. and Y. Kogan, 2000:  A New Cloud Physics
    # Parameterization in a Large-Eddy Simulation Model of Marine Stratocumulus.
    # Mon. Wea. Rev., 128, 229--243.
    #-----------------------------------------------------------------------
    Vrr = -(0.012 * (1.e6 * KK_mean_vol_rad) - 0.2)
    VNr = -(0.007 * (1.e6 * KK_mean_vol_rad) - 0.1)
    if l_clip_positive_sed:
        Vrr = jnp.minimum(Vrr, 0.0)
        VNr = jnp.minimum(VNr, 0.0)
        above = (jnp.arange(nzt)[None, :] > cloud_top_level[:, None]) & (cloud_top_level[:, None] > 0)
        Vrr = jnp.where(above, 0.0, Vrr)
        VNr = jnp.where(above, 0.0, VNr)
    return Vrr.at[:, -1].set(0.0), VNr.at[:, -1].set(0.0)

def KK_microphys_output(ngrdcol, nzt, hydromet_dim, hm_metadata, Vrr, VNr, rrm_mc, Nrm_mc):
    hydromet_mc = jnp.zeros(rrm_mc.shape + (hydromet_dim,))
    hydromet_vel = jnp.zeros_like(hydromet_mc)
    hydromet_mc = hydromet_mc.at[..., hm_metadata.iirr].set(rrm_mc)
    hydromet_mc = hydromet_mc.at[..., hm_metadata.iiNr].set(Nrm_mc)
    hydromet_vel = hydromet_vel.at[..., hm_metadata.iirr].set(Vrr)
    hydromet_vel = hydromet_vel.at[..., hm_metadata.iiNr].set(VNr)
    return hydromet_mc, hydromet_vel
