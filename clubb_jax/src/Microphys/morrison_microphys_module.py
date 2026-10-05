"""CLUBB–Morrison interface, following morrison_microphys_module.F90.

Columns are batched; out arguments are returned in source order. The WRF core's inout fields and
output diagnostics are returned under their Fortran names. The source's default-real conversions
are retained (float32 in the debug build). Atmospheric validation currently covers LBA's
warm-cloud/rain regime only.
"""

import jax
import jax.numpy as jnp
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, Ls, grav, sec_per_day
from clubb_jax.src.CLUBB_core import model_flags
from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K, T_in_K2thlm
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.Morrison_microphys.module_mp_graupel import M2005MICRO_GRAUPEL


# -----------------------------------------------------------------------------
def morrison_microphys_driver(
    gr, ngrdcol, dt, nzt,                     # In
    hydromet_dim, hm_metadata,                # In
    l_latin_hypercube, thlm, wm_zt, p_in_Pa,  # In
    exner, rho, cloud_frac, w_std_dev,        # In
    dzq, rcm, Ncm, chi, rvm, hydromet,        # In
    saturation_formula,                       # In
    sample_weight,                            # In
    stats,                                    # InOut
):
    """Apply Morrison microphysics to all columns and update diagnostics.

    l_latin_hypercube selects point inputs and weighted sample diagnostics; sample_weight is
    dimensionless. Tendencies and velocities are returned with updated statistics, preserving the
    established source-derived order.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        dt: Model timestep [s]
        nzt: Number of thermodynamic levels.
        hydromet_dim: Number of precipitating hydrometeor fields.
        hm_metadata: Hydrometeor/PDF names, zero-based species indices and tolerances.
        l_latin_hypercube: Whether we're using latin hypercube sampling
        thlm: Liquid potential temperature [K]
        wm_zt: Mean vertical velocity on the thermo grid [m/s]
        p_in_Pa: Pressure [Pa]
        exner: Exner function [-]
        rho: Density on thermo. grid [kg/m^3]
        cloud_frac: Cloud fraction [-]
        w_std_dev: Standard deviation of vertical vel. [m/s]
        dzq: Change in altitude [m]
        rcm: Cloud water mixing ratio [kg/kg]
        Ncm: Grid mean value for cloud droplet conc. [#/kg]
        chi: The variable 's' from Mellor [kg/kg]
        rvm: Vapor water mixing ratio [kg/kg]
        hydromet: Hydrometeor species [units vary]
        saturation_formula: Choice of liquid/ice saturation formula.
        sample_weight: SILHS sample weight; unity for an ordinary call
        stats: Immutable statistics state; return its updated value.
    """

    # Wrapper for the Morrison microphysics.
    # Description:
    # Wrapper for the Morrison microphysics
    #
    # References:
    # None
    # -----------------------------------------------------------------------
    # Fortran internal procedures use host association; nonlocal assignment
    # carries the immutable JaxStats object. Python definitions precede calls.
    def update_microphys_stat(name, value):
        nonlocal stats
        output_name = name

        # Non-interactive samples contribute only to the LH process diagnostics.
        if (
            l_latin_hypercube
            and parameters.lh_microphys_type == parameters.lh_microphys_non_interactive
        ):
            if name not in ("rrm_auto", "rrm_accr", "rrm_evap", "Nrm_auto", "Nrm_evap"):
                return
            output_name = "lh_" + name

        # est_silhs_tndcy performs the source sub-timestep average after sampling.
        if l_latin_hypercube:
            value = value * sample_weight
        stats = stats.update(output_name, value)

    def update_microphys_stat_sfc(name, value):
        # Non-interactive SILHS does not contribute to ordinary surface statistics.
        nonlocal stats
        if (
            l_latin_hypercube
            and parameters.lh_microphys_type == parameters.lh_microphys_non_interactive
        ):
            return
        if l_latin_hypercube:
            value = value * sample_weight[:, 0]
        stats = stats.update(name, value)

    def print_morr_error_output():
        # Host diagnostic adaptation: emit batched state and process mappings
        # instead of Fortran's per-level formatted field dump.
        nonfinite = jnp.any(
            jnp.stack(
                [
                    jnp.any(~jnp.isfinite(value))
                    for value in (
                        *final.values(),
                        diagnostic["rcm_mc"],
                        diagnostic["rvm_mc"],
                        diagnostic["T_mc"],
                    )
                ]
            )
        )
        jax.lax.cond(
            nonfinite,
            lambda _: jax.debug.print(
                "non-finite detected in a Morrison microphysics tendency\n"
                "altitude={z}\nfinal={final}\ndiagnostic={diagnostic}",
                z=gr.zt,
                final=final,
                diagnostic=diagnostic,
                ordered=True,
            ),
            lambda _: None,
            operand=None,
        )

    # The version of the Morrison 2005 microphysics that is in SAM.
    iirr, iiNr = hm_metadata.iirr, hm_metadata.iiNr
    iiri, iiNi = hm_metadata.iiri, hm_metadata.iiNi
    iirs, iiNs = hm_metadata.iirs, hm_metadata.iiNs
    iirg, iiNg = hm_metadata.iirg, hm_metadata.iiNg
    # Language adaptation of Fortran REAL(...), including its default-real
    # rounding before the numerical core and conversion back to core_rknd.
    zero = jnp.zeros_like(rcm)

    # Determine temperature.
    T_in_K = jnp.asarray(thlm2T_in_K(thlm, exner, rcm), dtype=jnp.float32)
    if l_latin_hypercube:
        # A sample represents one point: no grid-box cloud fraction or subgrid
        # velocity spread is supplied. Clip its updraft for droplet activation.
        w_thresh = 0.1  # Source minimum sample vertical velocity [m/s]
        cloud_frac_in = jnp.zeros_like(T_in_K)
        wm_zt_r4 = jnp.asarray(jnp.maximum(wm_zt, w_thresh), dtype=jnp.float32)
        w_std_dev_r4 = jnp.zeros_like(T_in_K)
    else:
        cloud_frac_in = jnp.asarray(cloud_frac, dtype=jnp.float32)
        wm_zt_r4 = jnp.asarray(wm_zt, dtype=jnp.float32)
        w_std_dev_r4 = jnp.asarray(w_std_dev, dtype=jnp.float32)
    rcm_r4, rvm_r4, Ncm_r4 = (
        jnp.asarray(rcm, dtype=jnp.float32),
        jnp.asarray(rvm, dtype=jnp.float32),
        jnp.asarray(Ncm, dtype=jnp.float32),
    )
    hydromet_r4 = jnp.asarray(hydromet, dtype=jnp.float32)

    # Below the homogeneous-freezing threshold, optionally remove liquid
    # cloud and cloud number before the core call (source 236.15 K threshold).
    if model_flags.l_evaporate_cold_rcm:
        cold = T_in_K < 236.15
        rcm_r4 = jnp.where(cold, 0.0, rcm_r4)
        cloud_frac_in = jnp.where(cold, 0.0, cloud_frac_in)
        Ncm_r4 = jnp.where(cold, 0.0, Ncm_r4)

    # Unpack hydrometeor arrays; absent species have zero extent in metadata.
    rrm, Nrm = hydromet[..., iirr], hydromet[..., iiNr]
    rim, Nim, rsm, Nsm, rgm, Ngm = (zero,) * 6
    species = [("qr", iirr), ("nr", iiNr)]
    if parameters.l_ice_microphys:
        rim, Nim = hydromet[..., iiri], hydromet[..., iiNi]
        rsm, Nsm = hydromet[..., iirs], hydromet[..., iiNs]
        species += [("qi", iiri), ("ni", iiNi), ("qs", iirs), ("ns", iiNs)]
        if parameters.l_graupel:
            rgm, Ngm = hydromet[..., iirg], hydromet[..., iiNg]
            species += [("qg", iirg), ("ng", iiNg)]

    # Moist energy and total water before microphysics, for conservation residuals.
    hl_before = (
        Cp * jnp.asarray(T_in_K, dtype=jnp.float64)
        + grav * gr.zt
        - Lv * (rcm + rrm)
        - Ls * (rim + rsm + rgm)
    )
    qto_before = rvm + rcm + rrm + rim + rsm + rgm

    # Call the one-column Morrison microphysics core. The leading batch axis
    # replaces Fortran's loop over i; output-only arguments are returned by name.
    rim_r4, rsm_r4 = jnp.asarray(rim, dtype=jnp.float32), jnp.asarray(rsm, dtype=jnp.float32)
    rrm_r4 = jnp.asarray(rrm, dtype=jnp.float32)
    Nim_r4, Nsm_r4 = jnp.asarray(Nim, dtype=jnp.float32), jnp.asarray(Nsm, dtype=jnp.float32)
    Nrm_r4 = jnp.asarray(Nrm, dtype=jnp.float32)
    P_in_pa_r4, rho_r4 = jnp.asarray(p_in_Pa, dtype=jnp.float32), jnp.asarray(
        rho, dtype=jnp.float32
    )
    dzq_r4 = jnp.asarray(dzq, dtype=jnp.float32)
    rgm_r4, Ngm_r4 = jnp.asarray(rgm, dtype=jnp.float32), jnp.asarray(Ngm, dtype=jnp.float32)
    output = M2005MICRO_GRAUPEL(
        rcm_r4, rim_r4, rsm_r4, rrm_r4, Ncm_r4, Nim_r4, Nsm_r4, Nrm_r4,      # In
        T_in_K, rvm_r4, P_in_pa_r4, rho_r4, dzq_r4, wm_zt_r4, w_std_dev_r4,  # In
        jnp.asarray(dt, dtype=jnp.float32),                                  # In
        1, 1, 1, 1, 1, nzt,                                                  # In
        1, 1, 1, 1, 1, nzt,                                                  # In
        rgm_r4, Ngm_r4,                                                      # In
        cloud_frac_in,                                                       # In
    )
    # JAX return adaptation: associate source inout/output arrays with the
    # existing CLUBB state/statistics names, without changing their precision.
    final = {
        key: output[source]
        for key, source in dict(
            qc="QC3D",
            qi="QI3D",
            qs="QNI3D",
            qr="QR3D",
            nc="NC3D",
            ni="NI3D",
            ns="NS3D",
            nr="NR3D",
            qg="QG3D",
            ng="NG3D",
            T="T3D",
            qv="QV3D",
        ).items()
    }
    diagnostic = dict(output)
    for name, source in dict(
        qc_sten="QCSTEN",
        qr_sten="QRSTEN",
        qi_sten="QISTEN",
        qs_sten="QNISTEN",
        qg_sten="QGSTEN",
        rcm_mc="QC3DTEN",
        rvm_mc="QV3DTEN",
        T_mc="T3DTEN",
        rain_vel="FR",
        eff_rad_cloud="EFFC",
        eff_rad_ice="EFFI",
        eff_rad_snow="EFFS",
        eff_rad_rain="EFFR",
        eff_rad_graupel="EFFG",
    ).items():
        diagnostic[name] = output[source]
    diagnostic["Morr_precip_rate"] = output["PRECRT"]
    diagnostic["Morr_snow_rate"] = output["SNOWRT"]
    final = {name: jnp.asarray(value, dtype=jnp.float64) for name, value in final.items()}
    diagnostic = {name: jnp.asarray(value, dtype=jnp.float64) for name, value in diagnostic.items()}
    rcm_sten, rrm_sten = diagnostic["qc_sten"], diagnostic["qr_sten"]
    rim_sten, rsm_sten, rgm_sten = (
        diagnostic["qi_sten"],
        diagnostic["qs_sten"],
        diagnostic["qg_sten"],
    )
    if clubb_at_least_debug_level(2):
        print_morr_error_output()

    # Subtract the core sedimentation tendencies from the energy/water change
    # to diagnose conservation independently of fallout from the column.
    hl_after = (
        Cp * final["T"]
        + grav * gr.zt
        - Lv * (final["qc"] + final["qr"])
        - Ls * (final["qi"] + final["qs"] + final["qg"])
    )
    hl_on_Cp_residual = (
        hl_after
        - hl_before
        - dt * Lv * (rcm_sten + rrm_sten)
        - dt * Ls * (rim_sten + rsm_sten + rgm_sten)
    ) / Cp
    qto_after = sum(final[name] for name in ("qv", "qc", "qr", "qi", "qs", "qg"))
    qto_residual = (
        qto_after - qto_before - dt * (rcm_sten + rrm_sten + rim_sten + rsm_sten + rgm_sten)
    )

    # Pack hydrometeor arrays.
    for name, index in species:
        hydromet_r4 = hydromet_r4.at[..., index].set(jnp.asarray(final[name], dtype=jnp.float32))
    rcm_mc, rvm_mc = diagnostic["rcm_mc"], diagnostic["rvm_mc"]
    rrm_auto, rrm_accr, rrm_evap = diagnostic["PRC"], diagnostic["PRA"], diagnostic["PRE"]
    Nrm_auto, Nrm_evap = diagnostic["NPRC1"], diagnostic["NSUBR"]

    # Finite differences include the core clipping in returned CLUBB tendencies.
    # Convert temperature/cloud liquid back to liquid-water potential temperature.
    hydromet_mc = (jnp.asarray(hydromet_r4, dtype=jnp.float64) - hydromet) / dt
    Ncm_mc = (final["nc"] - Ncm) / dt
    thlm_mc = (T_in_K2thlm(final["T"], exner, final["qc"]) - thlm) / dt
    # Sedimentation is handled within the Morrison microphysics.
    hydromet_vel_zt = jnp.zeros_like(hydromet)
    hydromet_vel_zt = hydromet_vel_zt.at[..., iirr].set(-diagnostic["rain_vel"])

    # Integrate the snow sedimentation tendency to check its column conservation.
    # Preserve the source's level-order accumulation independently per column.
    rsm_sd_morr_int = jax.lax.fori_loop(
        0,
        nzt,
        lambda k, integral: integral + rho[:, k] * rsm_sten[:, k] * gr.dzt[:, k],
        jnp.zeros((ngrdcol,), dtype=rho.dtype),
    )
    update_microphys_stat_sfc("rs_sd_morr_int", rsm_sd_morr_int)
    if clubb_at_least_debug_level(1):
        jax.lax.cond(
            jnp.any(rsm_sd_morr_int > jnp.max(rsm_sten, axis=1)),
            lambda _: jax.debug.print(
                "Warning: rsm_sd_morr was not conservative! "
                "rsm_sd_morr_verical_integr = {value}",
                value=rsm_sd_morr_int,
                ordered=True,
            ),
            lambda _: None,
            operand=None,
        )

    # Process rates, conservation residuals, and sedimentation tendencies.
    update_microphys_stat("rrm_auto", rrm_auto)
    update_microphys_stat("rrm_accr", rrm_accr)
    update_microphys_stat("rrm_evap", rrm_evap)
    update_microphys_stat("Nrm_auto", Nrm_auto)
    update_microphys_stat("Nrm_evap", Nrm_evap)
    update_microphys_stat("hl_on_Cp_residual", hl_on_Cp_residual)
    update_microphys_stat("qto_residual", qto_residual)
    update_microphys_stat("rgm_sd_morr", rgm_sten)
    update_microphys_stat("rrm_sd_morr", rrm_sten)
    update_microphys_stat("rsm_sd_morr", rsm_sten)
    update_microphys_stat("rim_sd_mg_morr", rim_sten)
    update_microphys_stat("rcm_sd_mg_morr", rcm_sten)

    # Named core process diagnostics, retained in source output order.
    update_microphys_stat("PRC", diagnostic["PRC"])
    update_microphys_stat("PRA", diagnostic["PRA"])
    update_microphys_stat("PRE", diagnostic["PRE"])
    update_microphys_stat("PSMLT", diagnostic["PSMLT"])
    update_microphys_stat("EVPMS", diagnostic["EVPMS"])
    update_microphys_stat("PRACS", diagnostic["PRACS"])
    update_microphys_stat("EVPMG", diagnostic["EVPMG"])
    update_microphys_stat("PRACG", diagnostic["PRACG"])
    update_microphys_stat("PGMLT", diagnostic["PGMLT"])
    update_microphys_stat("MNUCCC", diagnostic["MNUCCC"])
    update_microphys_stat("PSACWS", diagnostic["PSACWS"])
    update_microphys_stat("PSACWI", diagnostic["PSACWI"])
    update_microphys_stat("QMULTS", diagnostic["QMULTS"])
    update_microphys_stat("QMULTG", diagnostic["QMULTG"])
    update_microphys_stat("PSACWG", diagnostic["PSACWG"])
    update_microphys_stat("PGSACW", diagnostic["PGSACW"])
    update_microphys_stat("PRD", diagnostic["PRD"])
    update_microphys_stat("PRCI", diagnostic["PRCI"])
    update_microphys_stat("PRAI", diagnostic["PRAI"])
    update_microphys_stat("QMULTR", diagnostic["QMULTR"])
    update_microphys_stat("QMULTRG", diagnostic["QMULTRG"])
    update_microphys_stat("MNUCCD", diagnostic["MNUCCD"])
    update_microphys_stat("PRACI", diagnostic["PRACI"])
    update_microphys_stat("PRACIS", diagnostic["PRACIS"])
    update_microphys_stat("EPRD", diagnostic["EPRD"])
    update_microphys_stat("MNUCCR", diagnostic["MNUCCR"])
    update_microphys_stat("PIACR", diagnostic["PIACR"])
    update_microphys_stat("PIACRS", diagnostic["PIACRS"])
    update_microphys_stat("PGRACS", diagnostic["PGRACS"])
    update_microphys_stat("PRDS", diagnostic["PRDS"])
    update_microphys_stat("EPRDS", diagnostic["EPRDS"])
    update_microphys_stat("PSACR", diagnostic["PSACR"])
    update_microphys_stat("PRDG", diagnostic["PRDG"])
    update_microphys_stat("EPRDG", diagnostic["EPRDG"])

    # Number sedimentation and number-process diagnostics.
    update_microphys_stat("NGSTEN", diagnostic["NGSTEN"])
    update_microphys_stat("NRSTEN", diagnostic["NRSTEN"])
    update_microphys_stat("NISTEN", diagnostic["NISTEN"])
    update_microphys_stat("NSSTEN", diagnostic["NSSTEN"])
    update_microphys_stat("NCSTEN", diagnostic["NCSTEN"])
    update_microphys_stat("NPRC1", diagnostic["NPRC1"])
    update_microphys_stat("NRAGG", diagnostic["NRAGG"])
    update_microphys_stat("NPRACG", diagnostic["NPRACG"])
    update_microphys_stat("NSUBR", diagnostic["NSUBR"])
    update_microphys_stat("NSMLTR", diagnostic["NSMLTR"])
    update_microphys_stat("NGMLTR", diagnostic["NGMLTR"])
    update_microphys_stat("NPRACS", diagnostic["NPRACS"])
    update_microphys_stat("NNUCCR", diagnostic["NNUCCR"])
    update_microphys_stat("NIACR", diagnostic["NIACR"])
    update_microphys_stat("NIACRS", diagnostic["NIACRS"])
    update_microphys_stat("NGRACS", diagnostic["NGRACS"])
    update_microphys_stat("NSMLTS", diagnostic["NSMLTS"])
    update_microphys_stat("NSAGG", diagnostic["NSAGG"])
    update_microphys_stat("NPRCI", diagnostic["NPRCI"])
    update_microphys_stat("NSCNG", diagnostic["NSCNG"])
    update_microphys_stat("NSUBS", diagnostic["NSUBS"])
    update_microphys_stat("PCC", diagnostic["PCC"])
    update_microphys_stat("NNUCCC", diagnostic["NNUCCC"])
    update_microphys_stat("NPSACWS", diagnostic["NPSACWS"])
    update_microphys_stat("NPRA", diagnostic["NPRA"])
    update_microphys_stat("NPRC", diagnostic["NPRC"])
    update_microphys_stat("NPSACWI", diagnostic["NPSACWI"])
    update_microphys_stat("NPSACWG", diagnostic["NPSACWG"])
    update_microphys_stat("NPRAI", diagnostic["NPRAI"])
    update_microphys_stat("NMULTS", diagnostic["NMULTS"])
    update_microphys_stat("NMULTG", diagnostic["NMULTG"])
    update_microphys_stat("NMULTR", diagnostic["NMULTR"])
    update_microphys_stat("NMULTRG", diagnostic["NMULTRG"])
    update_microphys_stat("NNUCCD", diagnostic["NNUCCD"])
    update_microphys_stat("NSUBI", diagnostic["NSUBI"])
    update_microphys_stat("NGMLTG", diagnostic["NGMLTG"])
    update_microphys_stat("NSUBG", diagnostic["NSUBG"])
    update_microphys_stat("NACT", diagnostic["NACT"])

    # Size/positivity corrections and instantaneous core mass/number fields.
    update_microphys_stat("SIZEFIX_NR", diagnostic["SIZEFIX_NR"])
    update_microphys_stat("SIZEFIX_NC", diagnostic["SIZEFIX_NC"])
    update_microphys_stat("SIZEFIX_NI", diagnostic["SIZEFIX_NI"])
    update_microphys_stat("SIZEFIX_NS", diagnostic["SIZEFIX_NS"])
    update_microphys_stat("SIZEFIX_NG", diagnostic["SIZEFIX_NG"])
    update_microphys_stat("NEGFIX_NR", diagnostic["NEGFIX_NR"])
    update_microphys_stat("NEGFIX_NC", diagnostic["NEGFIX_NC"])
    update_microphys_stat("NEGFIX_NI", diagnostic["NEGFIX_NI"])
    update_microphys_stat("NEGFIX_NS", diagnostic["NEGFIX_NS"])
    update_microphys_stat("NEGFIX_NG", diagnostic["NEGFIX_NG"])
    update_microphys_stat("NIM_MORR_CL", diagnostic["NIM_MORR_CL"])
    update_microphys_stat("QC_INST", diagnostic["QC_INST"])
    update_microphys_stat("QR_INST", diagnostic["QR_INST"])
    update_microphys_stat("QI_INST", diagnostic["QI_INST"])
    update_microphys_stat("QS_INST", diagnostic["QS_INST"])
    update_microphys_stat("QG_INST", diagnostic["QG_INST"])
    update_microphys_stat("NC_INST", diagnostic["NC_INST"])
    update_microphys_stat("NR_INST", diagnostic["NR_INST"])
    update_microphys_stat("NI_INST", diagnostic["NI_INST"])
    update_microphys_stat("NS_INST", diagnostic["NS_INST"])
    update_microphys_stat("NG_INST", diagnostic["NG_INST"])

    # Temperature tendency and effective radii.
    update_microphys_stat("T_in_K_mc", diagnostic["T_mc"])
    update_microphys_stat("eff_rad_cloud", diagnostic["eff_rad_cloud"])
    update_microphys_stat("eff_rad_ice", diagnostic["eff_rad_ice"])
    update_microphys_stat("eff_rad_snow", diagnostic["eff_rad_snow"])
    update_microphys_stat("eff_rad_rain", diagnostic["eff_rad_rain"])
    update_microphys_stat("eff_rad_graupel", diagnostic["eff_rad_graupel"])
    # Core fallout is accumulated over dt; convert after promotion to core_rknd.
    update_microphys_stat_sfc("precip_rate_sfc", diagnostic["Morr_precip_rate"] * sec_per_day / dt)
    update_microphys_stat_sfc("morr_snow_rate", diagnostic["Morr_snow_rate"] * sec_per_day / dt)
    rrm_auto_diag, rrm_accr_diag, rrm_evap_diag = rrm_auto, rrm_accr, rrm_evap
    Nrm_auto_diag, Nrm_evap_diag = Nrm_auto, Nrm_evap
    return (
        stats,
        hydromet_mc,
        hydromet_vel_zt,
        Ncm_mc,
        rcm_mc,
        rvm_mc,
        thlm_mc,
        rrm_auto_diag,
        rrm_accr_diag,
        rrm_evap_diag,
        Nrm_auto_diag,
        Nrm_evap_diag,
    )
