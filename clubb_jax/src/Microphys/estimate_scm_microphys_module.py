"""Sampled microphysics from estimate_scm_microphys_module.F90.

Output-only arguments become functional returns; sample and vertical loops use JAX scan/array
operations. The microphysics procedure is static at trace time.
"""

import jax
import jax.numpy as jnp
from clubb_jax.src.Microphys import parameters_microphys
from clubb_jax.src.SILHS.math_utilities import compute_sample_mean
from clubb_jax.src.SILHS.latin_hypercube_driver_module import (
    copy_X_nl_into_hydromet_all_pts,
)
from clubb_jax.src.SILHS.lh_microphys_var_covar_module import (
    lh_microphys_var_covar_driver_api,
)
from clubb_jax.src.Microphys.silhs_category_variance_module import (
    silhs_category_variance_driver,
)
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust
from clubb_jax.src.Microphys.advance_microphys_module import get_cloud_top_level


# -----------------------------------------------------------------------------
def est_silhs_tndcy(
    gr, ngrdcol, dt, nzt, nzm, num_samples,                        # In
    pdf_dim, hydromet_dim, hm_metadata,                            # In
    X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights,  # In
    pdf_params, precip_fracs, p_in_Pa, exner, rho,                 # In
    dzq, hydromet, rcm,                                            # In
    lh_rt_clipped, lh_thl_clipped,                                 # In
    lh_rc_clipped, lh_rv_clipped,                                  # In
    lh_Nc_clipped,                                                 # In
    l_lh_instant_var_covar_src,                                    # In
    saturation_formula,                                            # In
    stats,                                                         # InOut
    microphys_sub,                                                 # In
):
    """Estimate the tendency of a microphysics scheme via latin hypercube sampling

    Arguments:
        gr: Grid metadata and coordinate arrays
        ngrdcol: Number of model columns
        dt: Model timestep [s]
        nzt: Number of thermodynamic vertical levels
        nzm: Number of momentum vertical levels
        num_samples: Number of calls to microphysics
        pdf_dim: Number of variates
        hydromet_dim: Number of precipitating hydrometeor fields.
        hm_metadata: Hydrometeor/PDF variable index metadata
        X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
        X_mixt_comp_all_levs: Mixture component of each sample
        lh_sample_point_weights: Weight for cloud weighted sampling
        pdf_params: The PDF parameters
        precip_fracs: Precipitation fractions [-]
        p_in_Pa: Pressure [Pa]
        exner: Exner function [-]
        rho: Density on thermo. grid [kg/m^3]
        dzq: Difference in height per gridbox [m]
        hydromet: Hydrometeor species [units vary]
        rcm: Mean liquid water mixing ratio [kg/kg]
        lh_rt_clipped: rt generated from silhs sample points
        lh_thl_clipped: thl generated from silhs sample points
        lh_rc_clipped: rc generated from silhs sample points
        lh_rv_clipped: rv generated from silhs sample points
        lh_Nc_clipped: Nc generated from silhs sample points
        l_lh_instant_var_covar_src: Produce instantaneous var/covar tendencies [-]
        saturation_formula: Integer that stores the saturation formula to be used
        stats: JAX statistics state; updated state is returned
        microphys_sub: Static microphysics procedure (local KK or Morrison)
    """
    # ---- Begin Code ----

    # Unpack sampled vertical velocity and Mellor (1977) s = chi [kg/kg].
    w_all_points = X_nl_all_levs[..., hm_metadata.iiPDF_w]
    chi_all_points = X_nl_all_levs[..., hm_metadata.iiPDF_chi]
    hydromet_all_points, Ncn_all_points = copy_X_nl_into_hydromet_all_pts(
        nzt, pdf_dim, num_samples, ngrdcol, X_nl_all_levs,  # In
        hydromet_dim, hm_metadata, hydromet,                # In
    )

    # Sampled schemes do not use grid-box cloud fraction or velocity spread.
    # Finite zeros replace the source unused_var sentinel.
    cloud_frac_unused = jnp.zeros_like(rcm)
    w_std_dev_unused = jnp.zeros_like(rcm)
    stats_before = stats

    def sample_microphysics(stats, sample):
        # Call the microphysics scheme to obtain one sample-point tendency.
        result = microphys_sub(
            gr, ngrdcol, dt, nzt, hydromet_dim, hm_metadata,          # In
            True, lh_thl_clipped[:, sample],                         # In
            w_all_points[:, sample], p_in_Pa, exner, rho,           # In
            cloud_frac_unused, w_std_dev_unused, dzq,               # In
            lh_rc_clipped[:, sample], lh_Nc_clipped[:, sample],     # In
            chi_all_points[:, sample], lh_rv_clipped[:, sample],    # In
            hydromet_all_points[:, sample], saturation_formula,     # In
            lh_sample_point_weights[:, sample],                     # In
            stats,                                                   # InOut
        )
        return result[0], result[1:]

    stats, all_tendencies = jax.lax.scan(sample_microphysics, stats, jnp.arange(num_samples))

    # Mirror stats_update(sub_timestep_average=.true.) through the shared API.
    stats = stats.average_subtimesteps(stats_before, num_samples)
    (
        lh_hydromet_mc_all,
        lh_hydromet_vel_all,
        lh_Ncm_mc_all,
        lh_rcm_mc_all,
        lh_rvm_mc_all,
        lh_thlm_mc_all,
        rrm_auto_diag_all,
        rrm_accr_diag_all,
        rrm_evap_diag_all,
        Nrm_auto_diag_all,
        Nrm_evap_diag_all,
    ) = tuple(jnp.swapaxes(x, 0, 1) for x in all_tendencies)

    # Compute variance/covariance tendencies, if requested.
    if parameters_microphys.l_var_covar_src:
        tendencies = lh_microphys_var_covar_driver_api(
            nzt, num_samples, ngrdcol, dt, lh_sample_point_weights,   # In
            pdf_params, lh_rt_clipped, lh_thl_clipped, w_all_points,  # In
            lh_rcm_mc_all, lh_rvm_mc_all, lh_thlm_mc_all,             # In
            l_lh_instant_var_covar_src,                               # In
        )

        # Convert thermodynamic-grid tendencies to the momentum grid.
        lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, lh_wpthlp_mc, lh_rtpthlp_mc = tuple(
            zt2zm(nzm, nzt, ngrdcol, gr, x) for x in tendencies
        )

        # Statistical sampling of the LH variance/covariance tendencies.
        for name, x in zip(
            (
                "lh_rtp2_mc",
                "lh_thlp2_mc",
                "lh_wprtp_mc",
                "lh_wpthlp_mc",
                "lh_rtpthlp_mc",
            ),
            (lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, lh_wpthlp_mc, lh_rtpthlp_mc),
        ):
            stats = stats.update(name, x)
    else:
        lh_rtp2_mc = lh_thlp2_mc = lh_wprtp_mc = lh_wpthlp_mc = lh_rtpthlp_mc = jnp.zeros(
            (ngrdcol, nzm)
        )

    # Grid-box averages of the sampled tendencies, velocities and diagnostics.
    lh_hydromet_vel = compute_sample_mean(
        nzt, num_samples, ngrdcol,                                # In
        lh_sample_point_weights[..., None], lh_hydromet_vel_all,  # In
    )
    lh_hydromet_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,                               # In
        lh_sample_point_weights[..., None], lh_hydromet_mc_all,  # In
    )
    lh_Ncm_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_Ncm_mc_all,  # In
    )
    lh_rcm_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_rcm_mc_all,  # In
    )
    lh_rvm_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,               # In
        lh_sample_point_weights, lh_rvm_mc_all,  # In
    )
    lh_thlm_mc = compute_sample_mean(
        nzt, num_samples, ngrdcol,                # In
        lh_sample_point_weights, lh_thlm_mc_all,  # In
    )
    rrm_auto_diag_avg = compute_sample_mean(
        nzt, num_samples, ngrdcol,                   # In
        lh_sample_point_weights, rrm_auto_diag_all,  # In
    )
    rrm_accr_diag_avg = compute_sample_mean(
        nzt, num_samples, ngrdcol,                   # In
        lh_sample_point_weights, rrm_accr_diag_all,  # In
    )
    rrm_evap_diag_avg = compute_sample_mean(
        nzt, num_samples, ngrdcol,                   # In
        lh_sample_point_weights, rrm_evap_diag_all,  # In
    )
    Nrm_auto_diag_avg = compute_sample_mean(
        nzt, num_samples, ngrdcol,                   # In
        lh_sample_point_weights, Nrm_auto_diag_all,  # In
    )
    Nrm_evap_diag_avg = compute_sample_mean(
        nzt, num_samples, ngrdcol,                   # In
        lh_sample_point_weights, Nrm_evap_diag_all,  # In
    )

    # Adjust means if l_silhs_KK_convergence_adj_mean is enabled.
    if parameters_microphys.l_silhs_KK_convergence_adj_mean:
        (
            lh_Vrr,
            lh_VNr,
            rrm_mc,
            Nrm_mc,
            lh_rvm_mc,
            lh_rcm_mc,
            lh_thlm_mc,
            lh_rrm_src_adj,
            lh_Nrm_src_adj,
            lh_rrm_evap_adj,
            lh_Nrm_evap_adj,
        ) = adjust_KK_src_means(
            dt, nzt, ngrdcol, exner, rcm, hydromet[..., hm_metadata.iirr],                   # In
            hydromet[..., hm_metadata.iiNr], hydromet,                                       # In
            hydromet_dim, hm_metadata.iiri, rrm_auto_diag_avg, rrm_accr_diag_avg,            # In
            rrm_evap_diag_avg,                                                               # In
            Nrm_auto_diag_avg, Nrm_evap_diag_avg,                                            # In
            lh_hydromet_vel[..., hm_metadata.iirr], lh_hydromet_vel[..., hm_metadata.iiNr],  # InOut
        )
        lh_hydromet_vel = (
            lh_hydromet_vel.at[..., hm_metadata.iirr]
            .set(lh_Vrr)
            .at[..., hm_metadata.iiNr]
            .set(lh_VNr)
        )
        lh_hydromet_mc = (
            lh_hydromet_mc.at[..., hm_metadata.iirr]
            .set(rrm_mc)
            .at[..., hm_metadata.iiNr]
            .set(Nrm_mc)
        )
        for name, x in zip(
            ("lh_rrm_src_adj", "lh_Nrm_src_adj", "lh_rrm_evap_adj", "lh_Nrm_evap_adj"),
            (lh_rrm_src_adj, lh_Nrm_src_adj, lh_rrm_evap_adj, lh_Nrm_evap_adj),
        ):
            stats = stats.update(name, x)

    # Invoke the category variance sampler if requested in the statistics list.
    if stats.l_sample and stats.var_on_stats_list("silhs_var_cat_1"):
        stats = silhs_category_variance_driver(
            ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata,  # In
            X_nl_all_levs, X_mixt_comp_all_levs, lh_hydromet_mc_all,        # In
            lh_sample_point_weights, pdf_params, precip_fracs,              # In
            stats,                                                          # InOut
        )
    return (
        stats,
        lh_hydromet_mc,
        lh_hydromet_vel,
        lh_Ncm_mc,
        lh_rvm_mc,
        lh_rcm_mc,
        lh_thlm_mc,
        lh_rtp2_mc,
        lh_thlp2_mc,
        lh_wprtp_mc,
        lh_wpthlp_mc,
        lh_rtpthlp_mc,
    )


# -----------------------------------------------------------------------------
def adjust_KK_src_means(
    dt, nzt, ngrdcol, exner, rcm, rrm, Nrm, hydromet,  # In
    hydromet_dim, iiri, rrm_auto, rrm_accr, rrm_evap,  # In
    Nrm_auto, Nrm_evap,                                # In
    lh_Vrr, lh_VNr,                                    # InOut
):
    """Adjusts the means of microphysics terms for KK microphysics by calling the KK microphysics
    adjustment subroutine for every model column.

    Source reference: CLUBB ticket 558.

    Arguments:
        dt: Model timestep [s]
        nzt: Number of thermodynamic levels.
        ngrdcol: Number of grid columns.
        exner: Exner function [-]
        rcm: Mean liquid water mixing ratio [kg/kg]
        rrm: Rain water mixing ratio [kg/kg]
        Nrm: Rain drop concentration [num/kg]
        hydromet: Hydrometeor fields used for cloud-top detection [units vary].
        hydromet_dim: Number of precipitating hydrometeor fields.
        iiri: Zero-based ice mixing-ratio index; negative if absent.
        rrm_auto: Mean change in rain due to autoconversion [(kg/kg)/s]
        rrm_accr: Mean change in rain due to accretion [(kg/kg)/s]
        rrm_evap: Mean change in rain due to evap [(kg/kg)/s]
        Nrm_auto: Mean change in Nrm due to autoconversion [(num/kg)/s]
        Nrm_evap: Mean change in Nrm due to evaporation [(num/kg)/s]
        lh_Vrr: Mean sedimentation velocity of < r_r > [m/s]
        lh_VNr: Mean sedimentation velocity of < N_r > [m/s]
    """
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, adj_terms = KK_microphys_adjust(
        dt, exner, rcm, rrm, Nrm,  # In
        rrm_evap, rrm_auto,        # In
        rrm_accr, Nrm_evap,        # In
        Nrm_auto, True,            # In
        True,                      # In
    )

    # Clip positive rain mixing-ratio and number sedimentation velocities.
    lh_Vrr = lh_Vrr.at[:, :-1].set(jnp.minimum(lh_Vrr[:, :-1], 0.0))
    lh_VNr = lh_VNr.at[:, :-1].set(jnp.minimum(lh_VNr[:, :-1], 0.0))

    # Mean sedimentation above cloud top must be zero.
    cloud_top_level = get_cloud_top_level(
        nzt, ngrdcol, rcm, hydromet,  # In
        hydromet_dim, iiri,           # In
    )
    above_cloud = (
        (jnp.arange(nzt)[None, :] > cloud_top_level[:, None])
        & (jnp.arange(nzt)[None, :] < nzt - 1)
        & (jnp.arange(nzt)[None, :] >= 1)
    )
    lh_Vrr = jnp.where(above_cloud, 0.0, lh_Vrr)
    lh_VNr = jnp.where(above_cloud, 0.0, lh_VNr)

    # Set the upper boundary tendencies to zero.
    rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc = tuple(
        x.at[:, -1].set(0.0) for x in (rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc)
    )
    return (lh_Vrr, lh_VNr, rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc, *adj_terms)
