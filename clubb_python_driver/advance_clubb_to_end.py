"""CLUBB time-loop driver utilities.

This module contains the timestep advancement logic that was split out from
the previous monolithic Python driver file.
"""

from numpy import asfortranarray as f_arr

import clubb_f2py
from clubb_python import clubb_api
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx
from clubb_python_driver.advance_clubb_core import advance_clubb_core as _advance_clubb_core_py


def advance_clubb_to_end(state: dict, l_stdout: bool = True):
    """Run the CLUBB time loop."""

    dt_main = state['dt_main']
    dt_rad = state['dt_rad']
    time_initial = state['time_initial']
    ifinal = state['ifinal']
    l_stats = state['l_stats']

    rad_interval = int(dt_rad / dt_main)

    for itime_idx in range(ifinal):
        itime = itime_idx + 1
        time_current = time_initial + (itime - 1) * dt_main
        l_sample = False
        l_last_sample = False

        # ── Stats: begin timestep ───────────────────────────────────────
        if l_stats:
            clubb_api.stats_begin_timestep(itime_idx)
            stats_cfg = clubb_api.get_stats_config()
            l_sample = bool(stats_cfg[7])
            l_last_sample = bool(stats_cfg[8])

        # ── Compute thvm ────────────────────────────────────────────────
        _calculate_thvm(state)

        # ── Prescribe forcings (case-specific) ──────────────────────────
        _prescribe_forcings(state, itime, l_sample=l_sample)

        # ── Add microphysical/radiative tendencies to forcings ─────────
        state['rtm_forcing'] = state['rtm_forcing'] + state['rcm_mc']
        state['thlm_forcing'] = state['thlm_forcing'] + state['thlm_mc'] + state['radht']

        # ── Radiation contribution to thlp2 ─────────────────────────────
        _calculate_thlp2_rad(state)

        # ── Advance CLUBB core ──────────────────────────────────────────
        state['l_sample'] = l_sample
        _advance_clubb_core(state)

        # ── Radiation ───────────────────────────────────────────────────
        l_rad_itime = (itime % rad_interval == 0) or (itime == 1)

        if l_rad_itime:
            _advance_radiation(
                state=state,
                time_current=time_current,
                l_sample=(l_stats and l_sample),
            )

        # ── Driver-owned stats updates (mirrors Fortran driver) ────────
        if l_stats and l_sample:
            state['Ncm'] = state['Nc_in_cloud'] * state['cloud_frac']
            clubb_api.stats_update("Ncm", state['Ncm'])
            clubb_api.stats_update("Nc_in_cloud", state['Nc_in_cloud'])

        # ── Stats: end timestep ─────────────────────────────────────────
        if l_stats and l_last_sample:
            # Match the Fortran driver: write the stats record at the model
            # time reached at the end of the current timestep, not one whole
            # output window later.
            stats_time = float(time_current + dt_main)
            state['err_info'] = clubb_api.stats_end_timestep(
                stats_time,
                err_info=state['err_info'],
            )

        # ── Update time ─────────────────────────────────────────────────
        time_current = time_initial + itime * dt_main

        if l_stdout:
            print(f"iteration: {itime:8d} / {ifinal:8d}"
                  f" -- time = {time_current:10.1f} / {state['time_final']:10.1f}")


def _calculate_thvm(state: dict):
    """Update virtual potential temperature diagnostic."""
    state['thvm'] = clubb_api.calculate_thvm(
        nzt=state['nzt'],
        ngrdcol=state['ngrdcol'],
        thlm=state['thlm'],
        rtm=state['rtm'],
        rcm=state['rcm'],
        exner=state['exner'],
        thv_ds_zt=state['thv_ds_zt'],
    )


def _calculate_thlp2_rad(state: dict):
    """Apply radiation contribution to thlp2 forcing when enabled."""
    if not state['l_calc_thlp2_rad']:
        return

    state['thlp2_forcing'] = clubb_api.calculate_thlp2_rad(
        gr=state['gr'],
        ngrdcol=state['ngrdcol'],
        nzm=state['nzm'],
        nzt=state['nzt'],
        rcm=state['rcm'],
        thlprcp=state['thlprcp'],
        radht=state['radht'],
        clubb_params=state['clubb_params'],
        thlp2_forcing=state['thlp2_forcing'],
    )


def _advance_clubb_core(state: dict):
    """Advance the CLUBB core using either the API wrapper or Python port."""
    # Default path: compiled Fortran core through the Python API adapter.
    #_advance_clubb_core_api(state)
    # To switch back to the translated Python port instead, replace the line
    # above with: _advance_clubb_core_python(state)
    _advance_clubb_core_python(state)


def _advance_clubb_core_python(state: dict):
    """Advance the CLUBB core using the translated Python port."""
    (
        state['um'],
        state['vm'],
        state['up3'],
        state['vp3'],
        state['thlm'],
        state['rtm'],
        state['rtp3'],
        state['thlp3'],
        state['wp3'],
        state['upwp'],
        state['vpwp'],
        state['up2'],
        state['vp2'],
        state['wprtp'],
        state['wpthlp'],
        state['rtp2'],
        state['thlp2'],
        state['rtpthlp'],
        state['wp2'],
        state['sclrm'],
        state['sclrp3'],
        state['wpsclrp'],
        state['sclrp2'],
        state['sclrprtp'],
        state['sclrpthlp'],
        state['p_in_Pa'],
        state['exner'],
        state['rcm'],
        state['cloud_frac'],
        state['wp2thvp'],
        state['wp2up'],
        state['wpthvp'],
        state['rtpthvp'],
        state['thlpthvp'],
        state['sclrpthvp'],
        state['wp2rtp'],
        state['wp2thlp'],
        state['wpup2'],
        state['wpvp2'],
        state['ice_supersat_frac'],
        state['uprcp'],
        state['vprcp'],
        state['rc_coef_zm'],
        state['wp4'],
        state['wp2up2'],
        state['wp2vp2'],
        state['um_pert'],
        state['vm_pert'],
        state['upwp_pert'],
        state['vpwp_pert'],
        state['edsclrm'],
        state['pdf_params'],
        state['pdf_params_zm'],
        state['pdf_implicit_coefs_terms'],
        state['err_info'],
        state['rcm_in_layer'],
        state['cloud_cover'],
        state['w_up_in_cloud'],
        state['w_down_in_cloud'],
        state['cloudy_updraft_frac'],
        state['cloudy_downdraft_frac'],
        state['wprcp_out'],
        state['invrs_tau_zm'],
        state['Kh_zt'],
        state['Kh_zm'],
        state['thlprcp'],
        state['Lscale'],
        state['_sigma_sqd_w'],
        state['_rc_coef'],
        state['_rcp2_zt'],
        state['_wprtp2'],
        state['_wpthlp2'],
        state['_wprtpthlp'],
        state['_wp2rcp'],
        state['_rtprcp'],
        state['_rcp2'],
        state['_skw_velocity'],
        state['_cloud_frac_zm'],
        state['_ice_supersat_frac_zm'],
        state['_rtm_zm'],
        state['_thlm_zm'],
        state['_rcm_zm'],
        state['_rcm_supersat_adj'],
        state['_wp2sclrp'],
        state['_wpsclrp2'],
        state['_sclrprcp'],
        state['_wpsclrprtp'],
        state['_wpsclrpthlp'],
    ) = _advance_clubb_core_py(
        gr=state['gr'],
        nzm=state['nzm'],
        nzt=state['nzt'],
        ngrdcol=state['ngrdcol'],
        l_implemented=False,
        dt=state['dt_main'],
        fcor=state['fcor'],
        fcor_y=state['fcor_y'],
        sfc_elevation=state['sfc_elevation'],
        hydromet_dim=state['hydromet_dim'],
        sclr_dim=state['sclr_dim'],
        sclr_tol=state['sclr_tol'],
        edsclr_dim=state['edsclr_dim'],
        sclr_idx=state['sclr_idx'],
        thlm_forcing=state['thlm_forcing'],
        rtm_forcing=state['rtm_forcing'],
        um_forcing=state['um_forcing'],
        vm_forcing=state['vm_forcing'],
        sclrm_forcing=state['sclrm_forcing'],
        edsclrm_forcing=state['edsclrm_forcing'],
        wprtp_forcing=state['wprtp_forcing'],
        wpthlp_forcing=state['wpthlp_forcing'],
        rtp2_forcing=state['rtp2_forcing'],
        thlp2_forcing=state['thlp2_forcing'],
        rtpthlp_forcing=state['rtpthlp_forcing'],
        wm_zm=state['wm_zm'],
        wm_zt=state['wm_zt'],
        wpthlp_sfc=state['wpthlp_sfc'],
        wprtp_sfc=state['wprtp_sfc'],
        upwp_sfc=state['upwp_sfc'],
        vpwp_sfc=state['vpwp_sfc'],
        p_sfc=state['p_sfc'],
        wpsclrp_sfc=state['wpsclrp_sfc'],
        wpedsclrp_sfc=state['wpedsclrp_sfc'],
        upwp_sfc_pert=state['upwp_sfc_pert'],
        vpwp_sfc_pert=state['vpwp_sfc_pert'],
        rtm_ref=state['rtm_ref'],
        thlm_ref=state['thlm_ref'],
        um_ref=state['um_ref'],
        vm_ref=state['vm_ref'],
        ug=state['ug'],
        vg=state['vg'],
        p_in_Pa=state['p_in_Pa'],
        rho_zm=state['rho_zm'],
        rho=state['rho'],
        exner=state['exner'],
        rho_ds_zm=state['rho_ds_zm'],
        rho_ds_zt=state['rho_ds_zt'],
        invrs_rho_ds_zm=state['invrs_rho_ds_zm'],
        invrs_rho_ds_zt=state['invrs_rho_ds_zt'],
        thv_ds_zm=state['thv_ds_zm'],
        thv_ds_zt=state['thv_ds_zt'],
        l_mix_rat_hm=state['l_mix_rat_hm'],
        rfrzm=state['rfrzm'],
        wphydrometp=state['wphydrometp'],
        wp2hmp=state['wp2hmp'],
        rtphmp_zt=state['rtphmp_zt'],
        thlphmp_zt=state['thlphmp_zt'],
        host_dx=state['host_dx'],
        host_dy=state['host_dy'],
        clubb_params=state['clubb_params'],
        nu_vert_res_dep=state['nu_vert_res_dep'],
        lmin=state['lmin'],
        mixt_frac_max_mag=state['mixt_frac_max_mag'],
        t0=state['T0'],
        ts_nudge=state['ts_nudge'],
        rtm_min=state['rtm_min'],
        rtm_nudge_max_altitude=state['rtm_nudge_max_altitude'],
        clubb_config_flags=state['flags'],
        l_sample=state.get('l_sample', False),
        um=state['um'],
        vm=state['vm'],
        upwp=state['upwp'],
        vpwp=state['vpwp'],
        up2=state['up2'],
        vp2=state['vp2'],
        up3=state['up3'],
        vp3=state['vp3'],
        thlm=state['thlm'],
        rtm=state['rtm'],
        wprtp=state['wprtp'],
        wpthlp=state['wpthlp'],
        wp2=state['wp2'],
        wp3=state['wp3'],
        rtp2=state['rtp2'],
        rtp3=state['rtp3'],
        thlp2=state['thlp2'],
        thlp3=state['thlp3'],
        rtpthlp=state['rtpthlp'],
        sclrm=state['sclrm'],
        sclrp2=state['sclrp2'],
        sclrp3=state['sclrp3'],
        sclrprtp=state['sclrprtp'],
        sclrpthlp=state['sclrpthlp'],
        wpsclrp=state['wpsclrp'],
        edsclrm=state['edsclrm'],
        rcm=state['rcm'],
        cloud_frac=state['cloud_frac'],
        wpthvp=state['wpthvp'],
        wp2thvp=state['wp2thvp'],
        wp2up=state['wp2up'],
        rtpthvp=state['rtpthvp'],
        thlpthvp=state['thlpthvp'],
        sclrpthvp=state['sclrpthvp'],
        wp2rtp=state['wp2rtp'],
        wp2thlp=state['wp2thlp'],
        uprcp=state['uprcp'],
        vprcp=state['vprcp'],
        rc_coef_zm=state['rc_coef_zm'],
        wp4=state['wp4'],
        wpup2=state['wpup2'],
        wpvp2=state['wpvp2'],
        wp2up2=state['wp2up2'],
        wp2vp2=state['wp2vp2'],
        ice_supersat_frac=state['ice_supersat_frac'],
        um_pert=state['um_pert'],
        vm_pert=state['vm_pert'],
        upwp_pert=state['upwp_pert'],
        vpwp_pert=state['vpwp_pert'],
        pdf_params=state['pdf_params'],
        pdf_params_zm=state['pdf_params_zm'],
        pdf_implicit_coefs_terms=state['pdf_implicit_coefs_terms'],
        err_info=state['err_info'],
    )


def _advance_clubb_core_api(state: dict):
    """Advance the CLUBB core through the compiled Fortran API wrapper."""
    (
        state['um'],
        state['vm'],
        state['up3'],
        state['vp3'],
        state['thlm'],
        state['rtm'],
        state['rtp3'],
        state['thlp3'],
        state['wp3'],
        state['upwp'],
        state['vpwp'],
        state['up2'],
        state['vp2'],
        state['wprtp'],
        state['wpthlp'],
        state['rtp2'],
        state['thlp2'],
        state['rtpthlp'],
        state['wp2'],
        state['sclrm'],
        state['sclrp3'],
        state['wpsclrp'],
        state['sclrp2'],
        state['sclrprtp'],
        state['sclrpthlp'],
        state['p_in_Pa'],
        state['exner'],
        state['rcm'],
        state['cloud_frac'],
        state['wp2thvp'],
        state['wp2up'],
        state['wpthvp'],
        state['rtpthvp'],
        state['thlpthvp'],
        state['sclrpthvp'],
        state['wp2rtp'],
        state['wp2thlp'],
        state['wpup2'],
        state['wpvp2'],
        state['ice_supersat_frac'],
        state['uprcp'],
        state['vprcp'],
        state['rc_coef_zm'],
        state['wp4'],
        state['wp2up2'],
        state['wp2vp2'],
        state['um_pert'],
        state['vm_pert'],
        state['upwp_pert'],
        state['vpwp_pert'],
        state['edsclrm'],
        state['pdf_params'],
        state['pdf_params_zm'],
        state['pdf_implicit_coefs_terms'],
        state['err_info'],
        state['rcm_in_layer'],
        state['cloud_cover'],
        state['w_up_in_cloud'],
        state['w_down_in_cloud'],
        state['cloudy_updraft_frac'],
        state['cloudy_downdraft_frac'],
        state['wprcp_out'],
        state['invrs_tau_zm'],
        state['Kh_zt'],
        state['Kh_zm'],
        state['thlprcp'],
        state['Lscale'],
    ) = clubb_api.advance_clubb_core(
        gr=state['gr'],
        nzm=state['nzm'],
        nzt=state['nzt'],
        ngrdcol=state['ngrdcol'],
        l_implemented=False,
        dt=state['dt_main'],
        fcor=state['fcor'],
        fcor_y=state['fcor_y'],
        sfc_elevation=state['sfc_elevation'],
        hydromet_dim=state['hydromet_dim'],
        sclr_dim=state['sclr_dim'],
        sclr_tol=state['sclr_tol'],
        edsclr_dim=state['edsclr_dim'],
        sclr_idx=state['sclr_idx'],
        thlm_forcing=state['thlm_forcing'],
        rtm_forcing=state['rtm_forcing'],
        um_forcing=state['um_forcing'],
        vm_forcing=state['vm_forcing'],
        sclrm_forcing=state['sclrm_forcing'],
        edsclrm_forcing=state['edsclrm_forcing'],
        wprtp_forcing=state['wprtp_forcing'],
        wpthlp_forcing=state['wpthlp_forcing'],
        rtp2_forcing=state['rtp2_forcing'],
        thlp2_forcing=state['thlp2_forcing'],
        rtpthlp_forcing=state['rtpthlp_forcing'],
        wm_zm=state['wm_zm'],
        wm_zt=state['wm_zt'],
        wpthlp_sfc=state['wpthlp_sfc'],
        wprtp_sfc=state['wprtp_sfc'],
        upwp_sfc=state['upwp_sfc'],
        vpwp_sfc=state['vpwp_sfc'],
        p_sfc=state['p_sfc'],
        wpsclrp_sfc=state['wpsclrp_sfc'],
        wpedsclrp_sfc=state['wpedsclrp_sfc'],
        upwp_sfc_pert=state['upwp_sfc_pert'],
        vpwp_sfc_pert=state['vpwp_sfc_pert'],
        rtm_ref=state['rtm_ref'],
        thlm_ref=state['thlm_ref'],
        um_ref=state['um_ref'],
        vm_ref=state['vm_ref'],
        ug=state['ug'],
        vg=state['vg'],
        p_in_Pa=state['p_in_Pa'],
        rho_zm=state['rho_zm'],
        rho=state['rho'],
        exner=state['exner'],
        rho_ds_zm=state['rho_ds_zm'],
        rho_ds_zt=state['rho_ds_zt'],
        invrs_rho_ds_zm=state['invrs_rho_ds_zm'],
        invrs_rho_ds_zt=state['invrs_rho_ds_zt'],
        thv_ds_zm=state['thv_ds_zm'],
        thv_ds_zt=state['thv_ds_zt'],
        l_mix_rat_hm=state['l_mix_rat_hm'],
        rfrzm=state['rfrzm'],
        wphydrometp=state['wphydrometp'],
        wp2hmp=state['wp2hmp'],
        rtphmp_zt=state['rtphmp_zt'],
        thlphmp_zt=state['thlphmp_zt'],
        host_dx=state['host_dx'],
        host_dy=state['host_dy'],
        clubb_params=state['clubb_params'],
        nu_vert_res_dep=state['nu_vert_res_dep'],
        lmin=state['lmin'],
        mixt_frac_max_mag=state['mixt_frac_max_mag'],
        t0=state['T0'],
        ts_nudge=state['ts_nudge'],
        rtm_min=state['rtm_min'],
        rtm_nudge_max_altitude=state['rtm_nudge_max_altitude'],
        clubb_config_flags=state['flags'],
        um=state['um'],
        vm=state['vm'],
        upwp=state['upwp'],
        vpwp=state['vpwp'],
        up2=state['up2'],
        vp2=state['vp2'],
        up3=state['up3'],
        vp3=state['vp3'],
        thlm=state['thlm'],
        rtm=state['rtm'],
        wprtp=state['wprtp'],
        wpthlp=state['wpthlp'],
        wp2=state['wp2'],
        wp3=state['wp3'],
        rtp2=state['rtp2'],
        rtp3=state['rtp3'],
        thlp2=state['thlp2'],
        thlp3=state['thlp3'],
        rtpthlp=state['rtpthlp'],
        sclrm=state['sclrm'],
        sclrp2=state['sclrp2'],
        sclrp3=state['sclrp3'],
        sclrprtp=state['sclrprtp'],
        sclrpthlp=state['sclrpthlp'],
        wpsclrp=state['wpsclrp'],
        edsclrm=state['edsclrm'],
        rcm=state['rcm'],
        cloud_frac=state['cloud_frac'],
        wpthvp=state['wpthvp'],
        wp2thvp=state['wp2thvp'],
        wp2up=state['wp2up'],
        rtpthvp=state['rtpthvp'],
        thlpthvp=state['thlpthvp'],
        sclrpthvp=state['sclrpthvp'],
        wp2rtp=state['wp2rtp'],
        wp2thlp=state['wp2thlp'],
        uprcp=state['uprcp'],
        vprcp=state['vprcp'],
        rc_coef_zm=state['rc_coef_zm'],
        wp4=state['wp4'],
        wpup2=state['wpup2'],
        wpvp2=state['wpvp2'],
        wp2up2=state['wp2up2'],
        wp2vp2=state['wp2vp2'],
        ice_supersat_frac=state['ice_supersat_frac'],
        um_pert=state['um_pert'],
        vm_pert=state['vm_pert'],
        upwp_pert=state['upwp_pert'],
        vpwp_pert=state['vpwp_pert'],
        pdf_params=state['pdf_params'],
        pdf_params_zm=state['pdf_params_zm'],
        pdf_implicit_coefs_terms=state['pdf_implicit_coefs_terms'],
        err_info=state['err_info'],
    )


def _prescribe_forcings(state: dict, itime: int, l_sample: bool = False):
    """Set forcings for the current timestep using Fortran prescribe_forcings."""
    del l_sample  # Sampling logic is handled internally by Fortran stats state.

    set_fortran_grid(state['gr'])
    set_fortran_sclr_idx(state['sclr_idx'])
    set_fortran_err_info(state['err_info'])
    result = clubb_f2py.f2py_prescribe_forcings(
        int(state['sclr_dim']), int(state['edsclr_dim']), str(state['runtype']), int(state['sfctype']),
        float(state['time_initial'] + (itime - 1) * state['dt_main']),
        float(state['time_initial']), float(state['dt_main']),
        f_arr(state['um']), f_arr(state['vm']), f_arr(state['thlm']), f_arr(state['p_in_Pa']),
        f_arr(state['exner']), f_arr(state['rho']), f_arr(state['rho_zm']), f_arr(state['thvm']),
        f_arr(state['gr'].zt),
        bool(state['l_t_dependent']), bool(state['l_ignore_forcings']), bool(state['l_input_xpwp_sfc']),
        bool(state['l_modify_bc_for_cnvg_test']),
        int(state['flags'].saturation_formula), bool(state['flags'].l_add_dycore_grid),
        int(state['flags'].grid_remap_method), int(state['flags'].grid_adapt_in_time_method),
        f_arr(state['rtm']), f_arr(state['wm_zm']), f_arr(state['wm_zt']), f_arr(state['ug']),
        f_arr(state['vg']), f_arr(state['um_ref']), f_arr(state['vm_ref']),
        f_arr(state['thlm_forcing']), f_arr(state['rtm_forcing']), f_arr(state['um_forcing']),
        f_arr(state['vm_forcing']), f_arr(state['wprtp_forcing']), f_arr(state['wpthlp_forcing']),
        f_arr(state['rtp2_forcing']), f_arr(state['thlp2_forcing']), f_arr(state['rtpthlp_forcing']),
        f_arr(state['wpsclrp']), f_arr(state['sclrm_forcing']), f_arr(state['edsclrm_forcing']),
        f_arr(state['wpthlp_sfc']), f_arr(state['wprtp_sfc']), f_arr(state['upwp_sfc']),
        f_arr(state['vpwp_sfc']), f_arr(state['T_sfc']), f_arr(state['p_sfc']),
        float(state['sens_ht']), float(state['latent_ht']),
        f_arr(state['wpsclrp_sfc']), f_arr(state['wpedsclrp_sfc']),
        nzm=int(state['nzm']), nzt=int(state['nzt']), ngrdcol=int(state['ngrdcol']),
    )

    (
        state['rtm'],
        state['wm_zm'],
        state['wm_zt'],
        state['ug'],
        state['vg'],
        state['um_ref'],
        state['vm_ref'],
        state['thlm_forcing'],
        state['rtm_forcing'],
        state['um_forcing'],
        state['vm_forcing'],
        state['wprtp_forcing'],
        state['wpthlp_forcing'],
        state['rtp2_forcing'],
        state['thlp2_forcing'],
        state['rtpthlp_forcing'],
        state['wpsclrp'],
        state['sclrm_forcing'],
        state['edsclrm_forcing'],
        state['wpthlp_sfc'],
        state['wprtp_sfc'],
        state['upwp_sfc'],
        state['vpwp_sfc'],
        state['T_sfc'],
        state['p_sfc'],
        state['sens_ht'],
        state['latent_ht'],
        state['wpsclrp_sfc'],
        state['wpedsclrp_sfc'],
        state['err_info'],
    ) = (*result, get_fortran_err_info())


def _advance_radiation(
    state: dict,
    time_current: float,
    l_sample: bool = False,
):
    """Advance radiation tendencies for currently supported schemes."""
    from clubb_python_driver.radiation import advance_radiation

    advance_radiation(state=state, time_current=time_current, l_sample=l_sample)
