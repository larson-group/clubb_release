"""CLUBB restart input from src/Input_fields/input_fields.F90.

Only the CLUBB statistics branch needed by restart_clubb is ported. General
LES inputfields remains gated in initialization. NetCDF reads/interpolation
are host I/O; restored arrays and immutable PDF parameters are returned.
JAX adaptation: read each stored column instead of duplicating the source
reader's first column. A single stored column can initialize several columns.
"""

from pathlib import Path

import numpy as np
import jax.numpy as jnp
from netCDF4 import Dataset, num2date

from clubb_jax.src.CLUBB_core.constants_clubb import eps, w_tol_sqd
from clubb_jax.src.CLUBB_core.interpolation import lin_interpolate_two_points


# Source inputfields module state, used only during host initialization.
stat_files = ()
clubb_day = 1
clubb_month = 1
clubb_year = 2000
# JAX radiation configuration is case-owned; initialization supplies the source
# soil/vegetation flag at this host-only inputfields boundary.
l_soil_veg = False
l_input_um = False
l_input_vm = False
l_input_rtm = False
l_input_thlm = False
l_input_wp2 = False
l_input_wprtp = False
l_input_wpthlp = False
l_input_wp3 = False
l_input_rtp2 = False
l_input_rtp3 = False
l_input_thlp2 = False
l_input_thlp3 = False
l_input_rtpthlp = False
l_input_upwp = False
l_input_vpwp = False
l_input_ug = False
l_input_vg = False
l_input_rcm = False
l_input_wm_zt = False
l_input_exner = False
l_input_em = False
l_input_p = False
l_input_rho = False
l_input_rho_zm = False
l_input_rho_ds_zm = False
l_input_rho_ds_zt = False
l_input_thv_ds_zm = False
l_input_thv_ds_zt = False
l_input_Lscale = False
l_input_Lscale_up = False
l_input_Lscale_down = False
l_input_Kh_zt = False
l_input_Kh_zm = False
l_input_tau_zm = False
l_input_tau_zt = False
l_input_wpthvp = False
l_input_wp2thvp = False
l_input_wp2up = False
l_input_rtpthvp = False
l_input_thlpthvp = False
l_input_wp2rtp = False
l_input_wp2thlp = False
l_input_uprcp = False
l_input_vprcp = False
l_input_rc_coef_zm = False
l_input_wp4 = False
l_input_wpup2 = False
l_input_wpvp2 = False
l_input_wp2up2 = False
l_input_wp2vp2 = False
l_input_iss_frac = False
l_input_radht = False
l_input_w_1 = False
l_input_w_2 = False
l_input_varnce_w_1 = False
l_input_varnce_w_2 = False
l_input_rt_1 = False
l_input_rt_2 = False
l_input_varnce_rt_1 = False
l_input_varnce_rt_2 = False
l_input_thl_1 = False
l_input_thl_2 = False
l_input_varnce_thl_1 = False
l_input_varnce_thl_2 = False
l_input_mixt_frac = False
l_input_chi_1 = False
l_input_chi_2 = False
l_input_stdev_chi_1 = False
l_input_stdev_chi_2 = False
l_input_rc_1 = False
l_input_rc_2 = False
l_input_w_1_zm = False
l_input_w_2_zm = False
l_input_varnce_w_1_zm = False
l_input_varnce_w_2_zm = False
l_input_mixt_frac_zm = False
l_input_thvm = False
l_input_rrm = False
l_input_Nrm = False
l_input_Ncm = False
l_input_rsm = False
l_input_rim = False
l_input_Nsm = False
l_input_Ngm = False
l_input_rgm = False
l_input_Nccnm = False
l_input_Nim = False
l_input_rrp2 = False
l_input_Nrp2 = False
l_input_wprrp = False
l_input_wpNrp = False
l_input_thlm_forcing = False
l_input_rtm_forcing = False
l_input_up2 = False
l_input_vp2 = False
l_input_sigma_sqd_w = False
l_input_cloud_frac = False
l_input_sigma_sqd_w_zt = False
l_input_veg_T_in_K = False
l_input_deep_soil_T_in_K = False
l_input_sfc_soil_T_in_K = False
l_input_wprtp_forcing = False
l_input_wpthlp_forcing = False
l_input_rtp2_forcing = False
l_input_thlp2_forcing = False
l_input_rtpthlp_forcing = False
l_input_thlprcp = False
l_input_rcm_mc = False
l_input_rvm_mc = False
l_input_thlm_mc = False
l_input_wprtp_mc = False
l_input_wpthlp_mc = False
l_input_rtp2_mc = False
l_input_thlp2_mc = False
l_input_rtpthlp_mc = False


# -----------------------------------------------------------------------------
def set_filenames(file_prefix):
    """Set the names of the netCDF files to be used for CLUBB restarts."""
    global stat_files

    stat_files = tuple(Path(str(file_prefix) + suffix) for suffix in (
        "_zt.nc", "_zm.nc", "_sfc.nc",
    ))
    if not stat_files[0].is_file():
        stat_files = (Path(str(file_prefix) + "_stats.nc"),) * 3
    return stat_files


# -----------------------------------------------------------------------------
def stat_fields_reader(
    gr, timestep, hydromet_dim, hm_metadata,                                               # In
    microphys_scheme, l_predict_Nc,                                                        # In
    um, upwp, vm, vpwp,                                                                    # InOut
    up2, vp2, rtm,                                                                         # InOut
    wprtp, thlm, wpthlp,                                                                   # InOut
    rtp2, rtp3,                                                                            # InOut
    thlp2, thlp3, rtpthlp,                                                                 # InOut
    wp2, wp3,                                                                              # InOut
    p_in_Pa, exner, rcm, cloud_frac,                                                       # InOut
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,                                             # InOut
    wp2rtp, wp2thlp, uprcp, vprcp,                                                         # InOut
    rc_coef_zm, wp4, wpup2,                                                                # InOut
    wpvp2, wp2up2,                                                                         # InOut
    wp2vp2, ice_supersat_frac,                                                             # InOut
    wm_zt, rho, rho_zm, rho_ds_zm,                                                         # InOut
    rho_ds_zt, thv_ds_zm, thv_ds_zt,                                                       # InOut
    thlm_forcing, rtm_forcing, wprtp_forcing,                                              # InOut
    wpthlp_forcing, rtp2_forcing,                                                          # InOut
    thlp2_forcing, rtpthlp_forcing,                                                        # InOut
    hydromet, hydrometp2, wphydrometp,                                                     # InOut
    Ncm, Nccnm, thvm, em,                                                                  # InOut
    tau_zm, tau_zt,                                                                        # InOut
    Kh_zt, Kh_zm, ug, vg,                                                                  # InOut
    thlprcp,                                                                               # InOut
    sigma_sqd_w, sigma_sqd_w_zt, radht,                                                    # InOut
    deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,                                         # InOut
    pdf_params, pdf_params_zm,                                                             # InOut
):
    """Read CLUBB statistics fields, in the source's thermo/momentum order.

    Arguments and returned inout fields retain the Fortran order. Profile
    arrays have shape (column, level); hydrometeors add a trailing species axis.

    Arguments (source order; profiles include all columns):
        gr: CLUBB thermodynamic/momentum grid [m].
        timestep: One-based output record to read [-].
        hydromet_dim: Number of hydrometeor species [-].
        hm_metadata: Hydrometeor/PDF species and variable indices [-].
        microphys_scheme: Selected microphysics scheme.
        l_predict_Nc: Whether cloud droplet concentration is prognostic [-].
        um: eastward grid-mean wind component (thermo. levs.)  [m/s]
        upwp: u'w' (momentum levels)                         [m^2/s^2]
        vm: northward grid-mean wind component (thermo. levs.) [m/s]
        vpwp: v'w' (momentum levels)                         [m^2/s^2]
        up2: u'^2 (momentum levels)                         [m^2/s^2]
        vp2: v'^2 (momentum levels)                         [m^2/s^2]
        rtm: total water mixing ratio, r_t (thermo. levels)     [kg/kg]
        wprtp: w' r_t' (momentum levels)                      [kg/kg m/s]
        thlm: liq. water pot. temp., th_l (thermo. levels)       [K]
        wpthlp: w'th_l' (momentum levels)                      [(m/s) K]
        rtp2: r_t'^2 (momentum levels)                       [(kg/kg)^2]
        rtp3: r_t'^3 (thermodynamic levels)                      [(kg/kg)^3]
        thlp2: th_l'^2 (momentum levels)                      [K^2]
        thlp3: th_l'^3 (thermodynamic levels)                     [K^3]
        rtpthlp: r_t'th_l' (momentum levels)                    [(kg/kg) K]
        wp2: w'^2 (momentum levels)                         [m^2/s^2]
        wp3: w'^3 (thermodynamic levels)                        [m^3/s^3]
        p_in_Pa: Air pressure (thermodynamic levels)                [Pa]
        exner: Exner function (thermodynamic levels)              [-]
        rcm: cloud water mixing ratio, r_c (thermo. levels)     [kg/kg]
        cloud_frac: cloud fraction (thermodynamic levels)              [-]
        wpthvp: < w' th_v' > (momentum levels)                 [kg/kg K]
        wp2thvp: < w'^2 th_v' > (thermodynamic levels)              [m^2/s^2 K]
        wp2up: < w'^2 u' > (thermodynamic levels)                 [m^3/s^3]
        rtpthvp: < r_t' th_v' > (momentum levels)               [kg/kg K]
        thlpthvp: < th_l' th_v' > (momentum levels)              [K^2]
        wp2rtp: w'^2 rt' (thermodynamic levels)                    [m^2/s^2 kg/kg]
        wp2thlp: w'^2 thl' (thermodynamic levels)                   [m^2/s^2 K]
        uprcp: < u' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
        vprcp: < v' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
        rc_coef_zm: Coef of X'r_c' in Eq. (34) (m-levs.)           [K/(kg/kg)]
        wp4: w'^4 (momentum levels)                         [m^4/s^4]
        wpup2: w'u'^2 (thermodynamic levels)                      [m^3/s^3]
        wpvp2: w'v'^2 (thermodynamic levels)                      [m^3/s^3]
        wp2up2: w'^2 u'^2 (momentum levels)                    [m^4/s^4]
        wp2vp2: w'^2 v'^2 (momentum levels)                    [m^4/s^4]
        ice_supersat_frac: ice cloud fraction (thermo. levels)                [-]
        wm_zt: vertical mean wind component on thermo. levels  [m/s]
        rho: Air density on thermodynamic levels             [kg/m^3]
        rho_zm: Air density on momentum levels                  [kg/m^3]
        rho_ds_zm: Dry, static density on momentum levels          [kg/m^3]
        rho_ds_zt: Dry, static density on thermo. levels           [kg/m^3]
        thv_ds_zm: Dry, base-state theta_v on momentum levels      [K]
        thv_ds_zt: Dry, base-state theta_v on thermo levels        [K]
        thlm_forcing: liquid potential temp. forcing (thermo. levels) [K/s]
        rtm_forcing: total water forcing (thermo. levels)      [(kg/kg)/s]
        wprtp_forcing: total water turbulent flux forcing (m-levs) [m*K/s^2]
        wpthlp_forcing: liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
        rtp2_forcing: total water variance forcing (m-levs)   [(kg/kg)^2/s]
        thlp2_forcing: liq pot temp variance forcing (m-levs)  [K^2/s]
        rtpthlp_forcing: <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]
        hydromet: Array of hydrometeors                [hm units]
        hydrometp2: Variance of a hydrometeor (m-levs.)  [<hm units>^2]
        wphydrometp: Covariance of w and a hydrometeor    [(m/s) <hm units>]
        Ncm: Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]
        Nccnm: Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
        thvm: Virtual potential temperature                        [K]
        em: Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
        tau_zm: Eddy dissipation time scale on momentum levels       [s]
        tau_zt: Eddy dissipation time scale on thermodynamic levels  [s]
        Kh_zt: Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
        Kh_zm: Eddy diffusivity coefficient on momentum levels      [m^2/s]
        ug: u geostrophic wind                                   [m/s]
        vg: v geostrophic wind                                   [m/s]
        thlprcp: thl'rc'                                              [K kg/kg]
        sigma_sqd_w: PDF width parameter (momentum levels)                [-]
        sigma_sqd_w_zt: PDF width parameter interpolated to t-levs.          [-]
        radht: SW + LW heating rate                                 [K/s]
        deep_soil_T_in_K: Deep soil temperature [K].
        sfc_soil_T_in_K: Surface soil temperature [K].
        veg_T_in_K: Vegetation temperature [K].
        pdf_params: PDF parameters (thermodynamic levels)    [units vary]
        pdf_params_zm: PDF parameters on momentum levels        [units vary]
    """
    iirr = hm_metadata.iirr
    iiNr = hm_metadata.iiNr
    iirs = hm_metadata.iirs
    iiri = hm_metadata.iiri
    iirg = hm_metadata.iirg
    iiNi = hm_metadata.iiNi
    # Normalize the mutable source hydrometeor array at the host I/O boundary.
    hydromet = jnp.asarray(hydromet)

    # -------------------------------------
    # CLUBB stats data
    # -------------------------------------
    # Thermo grid - zt file
    l_fatal_error = False
    tmp1 = jnp.zeros_like(um)

    um, l_read_error = get_clubb_variable_interpolated(
        l_input_um, stat_files[0], "um", gr.nzt, timestep,                                 # In
        gr.zt[0, :],                                                                       # In
        um,                                                                                # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    vm, l_read_error = get_clubb_variable_interpolated(
        l_input_vm, stat_files[0], "vm", gr.nzt, timestep,                                 # In
        gr.zt[0, :],                                                                       # In
        vm,                                                                                # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtm, l_read_error = get_clubb_variable_interpolated(
        l_input_rtm, stat_files[0], "rtm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        rtm,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlm, l_read_error = get_clubb_variable_interpolated(
        l_input_thlm, stat_files[0], "thlm", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        thlm,                                                                              # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp3, l_read_error = get_clubb_variable_interpolated(
        l_input_wp3, stat_files[0], "wp3", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        wp3,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tau_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_tau_zt, stat_files[0], "tau_zt", gr.nzt, timestep,                         # In
        gr.zt[0, :],                                                                       # In
        tau_zt,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_rrm, stat_files[0], "rrm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_rrm:
        hydromet = hydromet.at[:, :, iirr].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_rsm, stat_files[0], "rsm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_rsm:
        hydromet = hydromet.at[:, :, iirs].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_rim, stat_files[0], "rim", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_rim:
        hydromet = hydromet.at[:, :, iiri].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_rgm, stat_files[0], "rgm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_rgm:
        hydromet = hydromet.at[:, :, iirg].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    # Added variables for clubb_restart
    p_in_Pa, l_read_error = get_clubb_variable_interpolated(
        l_input_p, stat_files[0], "p_in_Pa", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        p_in_Pa,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    exner, l_read_error = get_clubb_variable_interpolated(
        l_input_exner, stat_files[0], "exner", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        exner,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    ug, l_read_error = get_clubb_variable_interpolated(
        l_input_ug, stat_files[0], "ug", gr.nzt, timestep,                                 # In
        gr.zt[0, :],                                                                       # In
        ug,                                                                                # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    vg, l_read_error = get_clubb_variable_interpolated(
        l_input_vg, stat_files[0], "vg", gr.nzt, timestep,                                 # In
        gr.zt[0, :],                                                                       # In
        vg,                                                                                # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rcm, l_read_error = get_clubb_variable_interpolated(
        l_input_rcm, stat_files[0], "rcm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        rcm,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wm_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_wm_zt, stat_files[0], "wm_zt", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        wm_zt,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rho, l_read_error = get_clubb_variable_interpolated(
        l_input_rho, stat_files[0], "rho", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        rho,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rho_ds_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_rho_ds_zt, stat_files[0], "rho_ds_zt", gr.nzt, timestep,                   # In
        gr.zt[0, :],                                                                       # In
        rho_ds_zt,                                                                         # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thv_ds_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_thv_ds_zt, stat_files[0], "thv_ds_zt", gr.nzt, timestep,                   # In
        gr.zt[0, :],                                                                       # In
        thv_ds_zt,                                                                         # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    Kh_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_Kh_zt, stat_files[0], "Kh_zt", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        Kh_zt,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thvm, l_read_error = get_clubb_variable_interpolated(
        l_input_thvm, stat_files[0], "thvm", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        thvm,                                                                              # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlm_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_thlm_forcing, stat_files[0], "thlm_forcing", gr.nzt, timestep,             # In
        gr.zt[0, :],                                                                       # In
        thlm_forcing,                                                                      # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtm_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_rtm_forcing, stat_files[0], "rtm_forcing", gr.nzt, timestep,               # In
        gr.zt[0, :],                                                                       # In
        rtm_forcing,                                                                       # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_Ncm, stat_files[0], "Ncm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_Ncm:
        Ncm = tmp1
    l_fatal_error = l_fatal_error or l_read_error

    Nccnm, l_read_error = get_clubb_variable_interpolated(
        l_input_Nccnm, stat_files[0], "Nccnm", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        Nccnm,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_Nim, stat_files[0], "Nim", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_Nim:
        hydromet = hydromet.at[:, :, iiNi].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    cloud_frac, l_read_error = get_clubb_variable_interpolated(
        l_input_cloud_frac, stat_files[0], "cloud_frac", gr.nzt, timestep,                 # In
        gr.zt[0, :],                                                                       # In
        cloud_frac,                                                                        # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_Nrm, stat_files[0], "Nrm", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        tmp1,                                                                              # InOut
    )
    if l_input_Nrm:
        hydromet = hydromet.at[:, :, iiNr].set(tmp1)
    l_fatal_error = l_fatal_error or l_read_error

    sigma_sqd_w_zt, l_read_error = get_clubb_variable_interpolated(
        l_input_sigma_sqd_w_zt, stat_files[0], "sigma_sqd_w_zt", gr.nzt, timestep,         # In
        gr.zt[0, :],                                                                       # In
        sigma_sqd_w_zt,                                                                    # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2thvp, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2thvp, stat_files[0], "wp2thvp", gr.nzt, timestep,                       # In
        gr.zt[0, :],                                                                       # In
        wp2thvp,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2up, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2up, stat_files[0], "wp2up", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        wp2up,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2rtp, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2rtp, stat_files[0], "wp2rtp", gr.nzt, timestep,                         # In
        gr.zt[0, :],                                                                       # In
        wp2rtp,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2thlp, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2thlp, stat_files[0], "wp2thlp", gr.nzt, timestep,                       # In
        gr.zt[0, :],                                                                       # In
        wp2thlp,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wpup2, l_read_error = get_clubb_variable_interpolated(
        l_input_wpup2, stat_files[0], "wpup2", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        wpup2,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wpvp2, l_read_error = get_clubb_variable_interpolated(
        l_input_wpvp2, stat_files[0], "wpvp2", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        wpvp2,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    ice_supersat_frac, l_read_error = get_clubb_variable_interpolated(
        l_input_iss_frac, stat_files[0], "ice_supersat_frac", gr.nzt, timestep,            # In
        gr.zt[0, :],                                                                       # In
        ice_supersat_frac,                                                                 # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    radht, l_read_error = get_clubb_variable_interpolated(
        l_input_radht, stat_files[0], "radht", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        radht,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    # PDF Parameters (needed for CLUBB restarts)
    # JAX adaptation: a profile scratch array replaces mutation of a PDF field.
    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_w_1, stat_files[0], "w_1", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.w_1,                                                                    # InOut
    )
    pdf_params = pdf_params.replace(w_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_w_2, stat_files[0], "w_2", gr.nzt, timestep,                               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.w_2,                                                                    # InOut
    )
    pdf_params = pdf_params.replace(w_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_w_1, stat_files[0], "varnce_w_1", gr.nzt, timestep,                 # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_w_1,                                                             # InOut
    )
    pdf_params = pdf_params.replace(varnce_w_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_w_2, stat_files[0], "varnce_w_2", gr.nzt, timestep,                 # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_w_2,                                                             # InOut
    )
    pdf_params = pdf_params.replace(varnce_w_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_rt_1, stat_files[0], "rt_1", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.rt_1,                                                                   # InOut
    )
    pdf_params = pdf_params.replace(rt_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_rt_2, stat_files[0], "rt_2", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.rt_2,                                                                   # InOut
    )
    pdf_params = pdf_params.replace(rt_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_rt_1, stat_files[0], "varnce_rt_1", gr.nzt, timestep,               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_rt_1,                                                            # InOut
    )
    pdf_params = pdf_params.replace(varnce_rt_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_rt_2, stat_files[0], "varnce_rt_2", gr.nzt, timestep,               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_rt_2,                                                            # InOut
    )
    pdf_params = pdf_params.replace(varnce_rt_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_thl_1, stat_files[0], "thl_1", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        pdf_params.thl_1,                                                                  # InOut
    )
    pdf_params = pdf_params.replace(thl_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_thl_2, stat_files[0], "thl_2", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        pdf_params.thl_2,                                                                  # InOut
    )
    pdf_params = pdf_params.replace(thl_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_thl_1, stat_files[0], "varnce_thl_1", gr.nzt, timestep,             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_thl_1,                                                           # InOut
    )
    pdf_params = pdf_params.replace(varnce_thl_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_thl_2, stat_files[0], "varnce_thl_2", gr.nzt, timestep,             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.varnce_thl_2,                                                           # InOut
    )
    pdf_params = pdf_params.replace(varnce_thl_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_mixt_frac, stat_files[0], "mixt_frac", gr.nzt, timestep,                   # In
        gr.zt[0, :],                                                                       # In
        pdf_params.mixt_frac,                                                              # InOut
    )
    pdf_params = pdf_params.replace(mixt_frac=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_chi_1, stat_files[0], "chi_1", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        pdf_params.chi_1,                                                                  # InOut
    )
    pdf_params = pdf_params.replace(chi_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_chi_2, stat_files[0], "chi_2", gr.nzt, timestep,                           # In
        gr.zt[0, :],                                                                       # In
        pdf_params.chi_2,                                                                  # InOut
    )
    pdf_params = pdf_params.replace(chi_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_stdev_chi_1, stat_files[0], "stdev_chi_1", gr.nzt, timestep,               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.stdev_chi_1,                                                            # InOut
    )
    pdf_params = pdf_params.replace(stdev_chi_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_stdev_chi_2, stat_files[0], "stdev_chi_2", gr.nzt, timestep,               # In
        gr.zt[0, :],                                                                       # In
        pdf_params.stdev_chi_2,                                                            # InOut
    )
    pdf_params = pdf_params.replace(stdev_chi_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_rc_1, stat_files[0], "rc_1", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.rc_1,                                                                   # InOut
    )
    pdf_params = pdf_params.replace(rc_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_rc_2, stat_files[0], "rc_2", gr.nzt, timestep,                             # In
        gr.zt[0, :],                                                                       # In
        pdf_params.rc_2,                                                                   # InOut
    )
    pdf_params = pdf_params.replace(rc_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    # Read in the zm file
    wp2, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2, stat_files[1], "wp2", gr.nzm, timestep,                               # In
        gr.zm[0, :],                                                                       # In
        wp2,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wprtp, l_read_error = get_clubb_variable_interpolated(
        l_input_wprtp, stat_files[1], "wprtp", gr.nzm, timestep,                           # In
        gr.zm[0, :],                                                                       # In
        wprtp,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wpthlp, l_read_error = get_clubb_variable_interpolated(
        l_input_wpthlp, stat_files[1], "wpthlp", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        wpthlp,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wpthvp, l_read_error = get_clubb_variable_interpolated(
        l_input_wpthvp, stat_files[1], "wpthvp", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        wpthvp,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtpthvp, l_read_error = get_clubb_variable_interpolated(
        l_input_rtpthvp, stat_files[1], "rtpthvp", gr.nzm, timestep,                       # In
        gr.zm[0, :],                                                                       # In
        rtpthvp,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlpthvp, l_read_error = get_clubb_variable_interpolated(
        l_input_thlpthvp, stat_files[1], "thlpthvp", gr.nzm, timestep,                     # In
        gr.zm[0, :],                                                                       # In
        thlpthvp,                                                                          # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtp2, l_read_error = get_clubb_variable_interpolated(
        l_input_rtp2, stat_files[1], "rtp2", gr.nzm, timestep,                             # In
        gr.zm[0, :],                                                                       # In
        rtp2,                                                                              # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlp2, l_read_error = get_clubb_variable_interpolated(
        l_input_thlp2, stat_files[1], "thlp2", gr.nzm, timestep,                           # In
        gr.zm[0, :],                                                                       # In
        thlp2,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtpthlp, l_read_error = get_clubb_variable_interpolated(
        l_input_rtpthlp, stat_files[1], "rtpthlp", gr.nzm, timestep,                       # In
        gr.zm[0, :],                                                                       # In
        rtpthlp,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    upwp, l_read_error = get_clubb_variable_interpolated(
        l_input_upwp, stat_files[1], "upwp", gr.nzm, timestep,                             # In
        gr.zm[0, :],                                                                       # In
        upwp,                                                                              # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    vpwp, l_read_error = get_clubb_variable_interpolated(
        l_input_vpwp, stat_files[1], "vpwp", gr.nzm, timestep,                             # In
        gr.zm[0, :],                                                                       # In
        vpwp,                                                                              # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    em, l_read_error = get_clubb_variable_interpolated(
        l_input_em, stat_files[1], "em", gr.nzm, timestep,                                 # In
        gr.zm[0, :],                                                                       # In
        em,                                                                                # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rho_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_rho_zm, stat_files[1], "rho_zm", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        rho_zm,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rho_ds_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_rho_ds_zm, stat_files[1], "rho_ds_zm", gr.nzm, timestep,                   # In
        gr.zm[0, :],                                                                       # In
        rho_ds_zm,                                                                         # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thv_ds_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_thv_ds_zm, stat_files[1], "thv_ds_zm", gr.nzm, timestep,                   # In
        gr.zm[0, :],                                                                       # In
        thv_ds_zm,                                                                         # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    Kh_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_Kh_zm, stat_files[1], "Kh_zm", gr.nzm, timestep,                           # In
        gr.zm[0, :],                                                                       # In
        Kh_zm,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    tau_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_tau_zm, stat_files[1], "tau_zm", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        tau_zm,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    up2, l_read_error = get_clubb_variable_interpolated(
        l_input_up2, stat_files[1], "up2", gr.nzm, timestep,                               # In
        gr.zm[0, :],                                                                       # In
        up2,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    vp2, l_read_error = get_clubb_variable_interpolated(
        l_input_vp2, stat_files[1], "vp2", gr.nzm, timestep,                               # In
        gr.zm[0, :],                                                                       # In
        vp2,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp4, l_read_error = get_clubb_variable_interpolated(
        l_input_wp4, stat_files[1], "wp4", gr.nzm, timestep,                               # In
        gr.zm[0, :],                                                                       # In
        wp4,                                                                               # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    uprcp, l_read_error = get_clubb_variable_interpolated(
        l_input_uprcp, stat_files[1], "uprcp", gr.nzm, timestep,                           # In
        gr.zm[0, :],                                                                       # In
        uprcp,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    vprcp, l_read_error = get_clubb_variable_interpolated(
        l_input_vprcp, stat_files[1], "vprcp", gr.nzm, timestep,                           # In
        gr.zm[0, :],                                                                       # In
        vprcp,                                                                             # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2up2, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2up2, stat_files[1], "wp2up2", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        wp2up2,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wp2vp2, l_read_error = get_clubb_variable_interpolated(
        l_input_wp2vp2, stat_files[1], "wp2vp2", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        wp2vp2,                                                                            # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    sigma_sqd_w, l_read_error = get_clubb_variable_interpolated(
        l_input_sigma_sqd_w, stat_files[1], "sigma_sqd_w", gr.nzm, timestep,               # In
        gr.zm[0, :],                                                                       # In
        sigma_sqd_w,                                                                       # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wprtp_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_wprtp_forcing, stat_files[1], "wprtp_forcing", gr.nzm, timestep,           # In
        gr.zm[0, :],                                                                       # In
        wprtp_forcing,                                                                     # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    wpthlp_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_wpthlp_forcing, stat_files[1], "wpthlp_forcing", gr.nzm, timestep,         # In
        gr.zm[0, :],                                                                       # In
        wpthlp_forcing,                                                                    # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtp2_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_rtp2_forcing, stat_files[1], "rtp2_forcing", gr.nzm, timestep,             # In
        gr.zm[0, :],                                                                       # In
        rtp2_forcing,                                                                      # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlp2_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_thlp2_forcing, stat_files[1], "thlp2_forcing", gr.nzm, timestep,           # In
        gr.zm[0, :],                                                                       # In
        thlp2_forcing,                                                                     # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rtpthlp_forcing, l_read_error = get_clubb_variable_interpolated(
        l_input_rtpthlp_forcing, stat_files[1], "rtpthlp_forcing", gr.nzm, timestep,       # In
        gr.zm[0, :],                                                                       # In
        rtpthlp_forcing,                                                                   # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    thlprcp, l_read_error = get_clubb_variable_interpolated(
        l_input_thlprcp, stat_files[1], "thlprcp", gr.nzm, timestep,                       # In
        gr.zm[0, :],                                                                       # In
        thlprcp,                                                                           # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    rc_coef_zm, l_read_error = get_clubb_variable_interpolated(
        l_input_rc_coef_zm, stat_files[1], "rc_coef_zm", gr.nzm, timestep,                 # In
        gr.zm[0, :],                                                                       # In
        rc_coef_zm,                                                                        # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    # PDF Parameters (needed for CLUBB restarts)
    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_w_1_zm, stat_files[1], "w_1_zm", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        pdf_params_zm.w_1,                                                                 # InOut
    )
    pdf_params_zm = pdf_params_zm.replace(w_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_w_2_zm, stat_files[1], "w_2_zm", gr.nzm, timestep,                         # In
        gr.zm[0, :],                                                                       # In
        pdf_params_zm.w_2,                                                                 # InOut
    )
    pdf_params_zm = pdf_params_zm.replace(w_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_w_1_zm, stat_files[1], "varnce_w_1_zm", gr.nzm, timestep,           # In
        gr.zm[0, :],                                                                       # In
        pdf_params_zm.varnce_w_1,                                                          # InOut
    )
    pdf_params_zm = pdf_params_zm.replace(varnce_w_1=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_varnce_w_2_zm, stat_files[1], "varnce_w_2_zm", gr.nzm, timestep,           # In
        gr.zm[0, :],                                                                       # In
        pdf_params_zm.varnce_w_2,                                                          # InOut
    )
    pdf_params_zm = pdf_params_zm.replace(varnce_w_2=profile)
    l_fatal_error = l_fatal_error or l_read_error

    profile, l_read_error = get_clubb_variable_interpolated(
        l_input_mixt_frac_zm, stat_files[1], "mixt_frac_zm", gr.nzm, timestep,             # In
        gr.zm[0, :],                                                                       # In
        pdf_params_zm.mixt_frac,                                                           # InOut
    )
    pdf_params_zm = pdf_params_zm.replace(mixt_frac=profile)
    l_fatal_error = l_fatal_error or l_read_error

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_veg_T_in_K, stat_files[2], "veg_T_in_K", 1, timestep,                      # In
        np.array([0.0]),                                                                   # In
        tmp1[:, :1],                                                                       # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    if l_input_veg_T_in_K:
        veg_T_in_K = tmp1[:, 0]

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_deep_soil_T_in_K, stat_files[2], "deep_soil_T_in_K", 1, timestep,          # In
        np.array([0.0]),                                                                   # In
        tmp1[:, :1],                                                                       # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    if l_input_deep_soil_T_in_K:
        deep_soil_T_in_K = tmp1[:, 0]

    tmp1, l_read_error = get_clubb_variable_interpolated(
        l_input_sfc_soil_T_in_K, stat_files[2], "sfc_soil_T_in_K", 1, timestep,            # In
        np.array([0.0]),                                                                   # In
        tmp1[:, :1],                                                                       # InOut
    )
    l_fatal_error = l_fatal_error or l_read_error

    if l_input_sfc_soil_T_in_K:
        sfc_soil_T_in_K = tmp1[:, 0]

    if l_fatal_error:
        raise ValueError("Failed to read CLUBB restart statistics")

    # Clipping on the variance of u, v and w.
    wp2 = jnp.maximum(wp2, w_tol_sqd)
    up2 = jnp.maximum(up2, w_tol_sqd)
    vp2 = jnp.maximum(vp2, w_tol_sqd)

    return (
        um, upwp, vm, vpwp,
        up2, vp2, rtm,
        wprtp, thlm, wpthlp,
        rtp2, rtp3,
        thlp2, thlp3, rtpthlp,
        wp2, wp3,
        p_in_Pa, exner, rcm, cloud_frac,
        wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
        wp2rtp, wp2thlp, uprcp, vprcp,
        rc_coef_zm, wp4, wpup2,
        wpvp2, wp2up2,
        wp2vp2, ice_supersat_frac,
        wm_zt, rho, rho_zm, rho_ds_zm,
        rho_ds_zt, thv_ds_zm, thv_ds_zt,
        thlm_forcing, rtm_forcing, wprtp_forcing,
        wpthlp_forcing, rtp2_forcing,
        thlp2_forcing, rtpthlp_forcing,
        hydromet, hydrometp2, wphydrometp,
        Ncm, Nccnm, thvm, em,
        tau_zm, tau_zt,
        Kh_zt, Kh_zm, ug, vg,
        thlprcp,
        sigma_sqd_w, sigma_sqd_w_zt, radht,
        deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
        pdf_params, pdf_params_zm
    )


# -----------------------------------------------------------------------------
def compute_timestep(filename, l_restart, time):
    """Determine the nearest saved output time; return a one-based record.

    Host adaptation of compute_timestep/open_netcdf_read: inspect actual time
    coordinates, including files whose output window starts after model start.
    Restart times must identify a saved record, rather than silently rounding.
    """
    with Dataset(filename) as ds:
        if "thlm" not in ds.variables:
            raise ValueError(f"Restart reference file has no thlm: {filename}")
        time_var = ds.variables["time"]
        dates = num2date(time_var[:], time_var.units,
                         calendar=getattr(time_var, "calendar", "standard"))
        # Model time is seconds since midnight on the configured start date.
        # Match open_netcdf_read's offset between file and model dates.
        model_midnight = dates[0].replace(
            year=clubb_year, month=clubb_month, day=clubb_day,
            hour=0, minute=0, second=0, microsecond=0,
        ) if len(dates) else None
        times = np.array([
            (date - model_midnight).total_seconds()
            for date in dates
        ])
        if not times.size or not np.isfinite(times).all():
            raise ValueError(f"Restart file has no valid output times: {filename}")
        nearest_timestep = int(np.argmin(np.abs(times - time)))
        if l_restart and abs(times[nearest_timestep] - time) > 1.0e-8:
            raise ValueError(f"time_restart={time} is not a saved output time in {filename}")
        return nearest_timestep + 1


# -----------------------------------------------------------------------------
def get_clubb_variable_interpolated(
    l_input_var, filename, varname, vardim, timestep,                                      # In
    clubb_heights,                                                                          # In
    variable_interpolated,                                                                  # InOut
):
    """Obtain a CLUBB profile and interpolate if needed (host I/O only).

    The source l_read_error output is returned with the restored array.
    NetCDF dimensions are normalized here, at the I/O boundary. Exact grids
    are copied directly to preserve restart bits.
    """
    if not l_input_var:
        return variable_interpolated, False

    try:
        with Dataset(filename) as ds:
            var = ds.variables[varname]
            # Preserve input_netcdf.get_var's source precision/read errors.
            if (
                np.dtype(var.dtype).kind != "f"
                or np.dtype(var.dtype).itemsize not in (4, 8)
            ):
                raise ValueError("Restart fields must use single or double precision")
            if "units" not in var.ncattrs() or not isinstance(var.units, str):
                raise ValueError(f"Restart field {varname} requires a units attribute")
            # CLUBB stats uses zero as _FillValue, including valid zero-valued
            # fields. Fortran reads raw values; disable netCDF4's masking too.
            var.set_auto_mask(False)
            var.set_auto_scale(False)
            selection = []
            grid_dim = None
            remaining_dims = []
            for dim in var.dimensions:
                if dim in ("time", "t", "T"):
                    if not 1 <= timestep <= len(ds.dimensions[dim]):
                        raise ValueError(f"Restart record {timestep} is outside {filename}")
                    selection.append(timestep - 1)
                elif dim in ("x", "X", "y", "Y", "lat", "lon", "latitude", "longitude"):
                    if len(ds.dimensions[dim]) != 1:
                        raise ValueError(f"Non-column horizontal dimension: {dim}")
                    selection.append(0)
                else:
                    selection.append(slice(None))
                    remaining_dims.append(dim)
                    if dim not in ("col", "column"):
                        grid_dim = dim
            data = var[tuple(selection)]
            data = np.asarray(data, dtype=np.float64)
            if grid_dim is None:
                data = data.reshape(-1, 1)
                heights = np.array([0.0])
            else:
                file_nz = data.shape[remaining_dims.index(grid_dim)]
                data = np.moveaxis(data, remaining_dims.index(grid_dim), -1)
                data = data.reshape(-1, file_nz)
                heights = (np.array([0.0]) if grid_dim == "sfc"
                           else np.asarray(ds.variables[grid_dim][:]))
            units = var.units.strip()
            if units == "g/kg":
                data = data / 1000.0
            elif units == "K/day":
                data = data / 86400.0
            elif units == "W/m2":
                raise ValueError("Cannot convert W/m2 restart field to MKS")

            if data.shape[0] not in (1, variable_interpolated.shape[0]):
                raise ValueError(f"Restart column count differs for {varname}")
            if heights.shape != (vardim,) or np.any(
                np.abs(heights - clubb_heights) > np.abs(heights + clubb_heights) * eps / 2
            ):
                # Source stat_file_average: linear interpolation in domain;
                # out-of-domain levels are infinite, except the lower ghost.
                order = np.argsort(heights)
                file_variable = data
                upper_lev_idx = np.searchsorted(heights[order], clubb_heights)
                upper_lev_idx = np.clip(upper_lev_idx, 1, max(len(heights) - 1, 1))
                lower_lev_idx = upper_lev_idx - 1
                if len(heights) > 1:
                    data = lin_interpolate_two_points(
                        clubb_heights, heights[order[upper_lev_idx]], heights[order[lower_lev_idx]],  # In
                        file_variable[:, order[upper_lev_idx]], file_variable[:, order[lower_lev_idx]],  # In
                    )
                else:
                    data = np.broadcast_to(file_variable, (file_variable.shape[0], vardim)).copy()
                # Exact levels retain their stored bits, as in the source.
                exact_lev_idx = np.clip(np.searchsorted(heights[order], clubb_heights), 0, len(heights) - 1)
                data = np.where(
                    clubb_heights == heights[order[exact_lev_idx]],
                    file_variable[:, order[exact_lev_idx]], data,
                )
                in_domain = (clubb_heights >= heights[order][0]) & (clubb_heights <= heights[order][-1])
                data = np.where(in_domain, data, np.inf)
                if clubb_heights[0] < 0.0 and clubb_heights[0] < heights[order][0]:
                    data[:, 0] = file_variable[:, order[0]]
            data = np.broadcast_to(data, variable_interpolated.shape)
            if not np.isfinite(data).all():
                raise ValueError(f"Restart grid does not cover {varname}")
            return jnp.asarray(data), False
    except (OSError, KeyError, ValueError, IndexError) as exc:
        print(f"Error reading {varname} from {filename}: {exc}")
        return variable_interpolated, True
