"""CLUBB time-loop driver utilities.

Mirrors the timestep loop of the Fortran `advance_clubb_to_end` subroutine (clubb_driver.F90): per step it
applies forcings, adds radiation tendencies, calls the closure advance (advance_clubb_core), advances
radiation, and accumulates stats. Split out from the previous monolithic Python driver file.
"""

import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.advance_clubb_core_module import advance_clubb_core
from clubb_jax.src.CLUBB_core.advance_helper_module import calculate_thlp2_rad
from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Benchmark_cases.prescribe_forcings import prescribe_forcings
from clubb_jax.src.Radiation.radiation_module import advance_clubb_radiation


def _err_code_summary(err_info) -> str:
    if err_info is None:
        return "err_info=None"
    err_code = getattr(err_info, "err_code", None)
    reason_code = getattr(err_info, "reason_code", None)
    reason_messages = err_info.reason_messages_host()
    pieces = []
    if err_code is None:
        pieces.append("err_code=None")
    else:
        pieces.append(f"err_code={err_code}")
    if reason_code is not None:
        pieces.append(f"reason_code={reason_code}")
    if reason_messages:
        pieces.append("message=" + " | ".join(reason_messages))
    return "; ".join(pieces)


def advance_clubb_to_end(state: dict, l_stdout: bool = True, max_steps: int | None = None):
    """Run the CLUBB time loop."""

    dt_main = state['dt_main']
    dt_rad = state['dt_rad']
    time_initial = state['time_initial']
    ifinal = state['ifinal']
    l_stats = state['l_stats']

    rad_interval = int(dt_rad / dt_main)
    n_steps = ifinal if max_steps is None else min(ifinal, max_steps)

    for itime_idx in range(n_steps):
        itime = itime_idx + 1
        time_current = time_initial + (itime - 1) * dt_main

        # ── Stats: begin timestep ───────────────────────────────────────
        # Begin stats collection and initialize JAX stats banks for core updates.
        if l_stats:
            stats_writer = state['stats_writer']
            l_sample, l_last_sample = stats_writer.begin_timestep(itime_idx)
            jax_stats = state.get('_jax_stats')
            if jax_stats is None:
                jax_stats = JaxStats.from_layout(
                    stats_writer.get_jax_layout(),
                    ncol=state['ngrdcol'],
                )
            state['_jax_stats'] = jax_stats.begin_timestep(
                l_sample=l_sample,
                reset_accumulators=stats_writer.l_reset,
            )
        else:
            state['_jax_stats'] = JaxStats.empty(
                l_sample=False,
                names=(),
                ncol=state['ngrdcol'],
                max_nlev=max(state['nzm'], state['nzt'], 1),
            )
            l_sample = False
            l_last_sample = False

        # ── Compute thvm ────────────────────────────────────────────────
        # Update virtual potential temperature diagnostic.
        state['thvm'] = calculate_thvm(
            nzt=state['nzt'],
            ngrdcol=state['ngrdcol'],
            thlm=state['thlm'],
            rtm=state['rtm'],
            rcm=state['rcm'],
            exner=state['exner'],
            thv_ds_zt=state['thv_ds_zt'],
        )

        # ── Prescribe forcings (case-specific) ──────────────────────────
        _prescribe_forcings(state, time_current)

        # ── Add radiation tendency to forcing ─────────────────────────────────
        state['thlm_forcing'] = state['thlm_forcing'] + state['radht']

        # ── Radiation contribution to thlp2 ─────────────────────────────
        # Apply radiation contribution to thlp2 forcing when enabled.
        if state['l_calc_thlp2_rad']:
            state['thlp2_forcing'] = calculate_thlp2_rad(
                state['ngrdcol'],
                state['nzm'],
                state['nzt'],
                state['gr'],
                state['rcm'],
                state['thlprcp'],
                state['radht'],
                state['clubb_params'],
                state['thlp2_forcing'],
            )

        # ── Advance CLUBB core ──────────────────────────────────────────
        state['l_sample'] = l_sample
        _advance_clubb_core(state)

        # ── Radiation ───────────────────────────────────────────────────
        l_rad_itime = (itime % rad_interval == 0) or (itime == 1)

        _advance_radiation(state=state, time_current=time_current, l_rad_itime=l_rad_itime)
        if state['err_info'].is_fatal():
            raise RuntimeError(
                "Fatal error in radiation; "
                f"{_err_code_summary(state['err_info'])}"
            )

        if l_stats and l_sample:
            state['Ncm'] = state['Nc_in_cloud'] * state['cloud_frac']
            state['_jax_stats'] = state['_jax_stats'].update("Ncm", state['Ncm'])
            state['_jax_stats'] = state['_jax_stats'].update(
                "Nc_in_cloud", state['Nc_in_cloud']
            )

        # ── Stats: end timestep ─────────────────────────────────────────
        if l_last_sample:
            # Finalize a sampled stats timestep through the F2PY-free writer.
            # ``time_current`` is the start of the model step.  The output timestamp
            # is the end of this step, not the end of the potentially multi-step
            # sampling interval (clubb_driver.F90, d4fbf466f).
            stats_time = float(time_current + state['dt_main'])
            state['_jax_stats'] = state['stats_writer'].end_timestep(
                stats_time,
                jax_stats=state['_jax_stats'],
            )

        # ── Update time ─────────────────────────────────────────────────
        time_current = time_initial + itime * dt_main
        if l_stdout:
            print(f"iteration: {itime:8d} / {ifinal:8d}"
                  f" -- time = {time_current:10.1f} / {state['time_final']:10.1f}")


def _prescribe_forcings(state: dict, time_current):
    """Apply case-specific forcings using the JAX translated routine."""
    (
        state['_jax_stats'],
        state['rtm'], state['wm_zm'], state['wm_zt'],
        state['ug'], state['vg'], state['um_ref'], state['vm_ref'],
        state['thlm_forcing'], state['rtm_forcing'], state['um_forcing'],
        state['vm_forcing'], state['wprtp_forcing'], state['wpthlp_forcing'],
        state['rtp2_forcing'], state['thlp2_forcing'], state['rtpthlp_forcing'],
        state['wpsclrp'], state['sclrm_forcing'], state['edsclrm_forcing'],
        state['wpthlp_sfc'], state['wprtp_sfc'],
        state['upwp_sfc'], state['vpwp_sfc'],
        state['T_sfc'], state['p_sfc'], state['sens_ht'], state['latent_ht'],
        state['wpsclrp_sfc'], state['wpedsclrp_sfc'], state['err_info'],
    ) = prescribe_forcings(
        state['gr'], state['nzm'], state['nzt'], state['ngrdcol'],
        state['sclr_dim'], state['edsclr_dim'], state['sclr_idx'],
        state['runtype'], state['sfctype'],
        jnp.asarray(time_current), state['time_initial'], state['dt_main'],
        state['um'], state['vm'], state['thlm'],
        state['p_in_Pa'], state['exner'], state['rho'],
        state['rho_zm'], state['thvm'],
        state['veg_T_in_K'],
        state['l_modify_bc_for_cnvg_test'],
        state['flags'].saturation_formula,
        state['_jax_stats'],
        state['rtm'], state['wm_zm'], state['wm_zt'],
        state['ug'], state['vg'], state['um_ref'], state['vm_ref'],
        state['thlm_forcing'], state['rtm_forcing'], state['um_forcing'],
        state['vm_forcing'], state['wprtp_forcing'], state['wpthlp_forcing'],
        state['rtp2_forcing'], state['thlp2_forcing'], state['rtpthlp_forcing'],
        state['wpsclrp'], state['sclrm_forcing'], state['edsclrm_forcing'],
        state['wpthlp_sfc'], state['wprtp_sfc'],
        state['upwp_sfc'], state['vpwp_sfc'],
        state['T_sfc'], state['p_sfc'], state['sens_ht'], state['latent_ht'],
        state['wpsclrp_sfc'], state['wpedsclrp_sfc'], state['err_info'],
    )


def _advance_clubb_core(state: dict):
    """Advance the CLUBB core using the JAX translated core."""
    core_kwargs = dict(
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
        stats=state['_jax_stats'],
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
    result = advance_clubb_core(**core_kwargs)
    if result is None:
        raise RuntimeError(
            "advance_clubb_core returned without outputs; "
            f"{_err_code_summary(state.get('err_info'))}"
        )

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
        state['_jax_stats'],
    ) = result
    if state['err_info'].is_fatal():
        nonfinite = []
        for name, value in state.items():
            if not hasattr(value, "dtype"):
                continue
            array = np.asarray(value)
            if array.dtype.kind in "fc" and not np.isfinite(array).all():
                nonfinite.append(
                    f"{name}(nan={int(np.isnan(array).sum())},inf={int(np.isinf(array).sum())})"
                )
        message = (
            "advance_clubb_core returned fatal err_info; "
            f"{_err_code_summary(state['err_info'])}"
        )
        if nonfinite:
            message += "; nonfinite=" + ", ".join(nonfinite)
        print(message)
        raise RuntimeError(message)


def _advance_radiation(
    state: dict,
    time_current: float,
    l_rad_itime: bool,
):
    """Advance radiation through the Fortran-named functional JAX driver.

    The source owns radiation arrays in ``radiation_variables_module``. The
    JAX port passes and unpacks those explicit functional values in source
    ownership order after the source argument list.
    """
    (
        state['_jax_stats'], state['err_info'],
        state['deep_soil_T_in_K'], state['sfc_soil_T_in_K'], state['veg_T_in_K'],
        state['radht'], state['Frad'],
        state['Frad_SW_up'], state['Frad_LW_up'],
        state['Frad_SW_down'], state['Frad_LW_down'],
        state['radht_SW'], state['radht_LW'], state['Frad_SW'], state['Frad_LW'],
    ) = advance_clubb_radiation(
        state['gr'], state['ngrdcol'], state['hydromet_dim'], 0, 0,
        jnp.asarray(l_rad_itime), state['dt_main'], state['day'], state['month'], state['year'],
        state['lat_vals'], state['lon_vals'], jnp.asarray(time_current), state['time_initial'],
        state['rho'], state['rho_zm'], state['p_in_Pa'], state['exner'],
        state['wpthlp_sfc'], state['wprtp_sfc'], state['p_sfc'],
        state['cloud_frac'], state['ice_supersat_frac'],
        state['thlm'], state['rtm'], state['rcm'],
        state['X_nl_all_levs'], state['lh_rt_clipped'], state['lh_thl_clipped'], state['lh_rc_clipped'],
        state['lh_sample_point_weights'], state['hydromet'], state['hm_metadata'],
        state['_jax_stats'], state['err_info'],
        state['deep_soil_T_in_K'], state['sfc_soil_T_in_K'], state['veg_T_in_K'],
        state['radht'], state['Frad'],
        state['Frad_SW_up'], state['Frad_LW_up'], state['Frad_SW_down'], state['Frad_LW_down'],
        state['radht_SW'], state['radht_LW'], state['Frad_SW'], state['Frad_LW'],
        state['radiation_parameters'],
    )
