"""JAX port of selected routines from ``src/CLUBB_core/numerical_check.F90``.

*_check subroutines were added to ensure that the
subroutines they are checking perform correctly
Joshua Fasching February 2008

rad_clipping has been replaced by rad_check as the new
subroutine only reports if there are invalid values.
Joshua Fasching March 2008

Porting deviations:
- The Fortran module reports diagnostics through formatted writes and mutates
  `err_info`; JAX helper predicates return booleans where the caller only needs
  validity and uses `ErrInfo` updates only for `parameterization_check`,
  `check_negative`, and `check_nan`.
- The Fortran generic `check_nan` has scalar, 1D, and 2D implementations. JAX
  uses one broadcasting implementation because `jnp.isfinite` handles all
  ranks used by the port.
- `check_clubb_settings_api` ports the parameter checks and the
  l_do_expldiff_rtm_thlm eddy-scalar checks. Its other flag checks are left
  out because the JAX driver rejects unsupported runtime options elsewhere.
  The model flag constants used by the Fortran routine are still mirrored
  below so dispatch values remain source-grounded.

Observed Fortran differences retained in this diagnostic pass:
- parameterization_check keeps evaluating after fatal NaNs instead of taking the Fortran early return; cold-level output is suppressed after NaNs.
- Full mixed-flag configuration checks and the detailed PDF validator schema remain unported; reduced host validity predicates do not carry error-header metadata.
"""


import sys

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp
import jax
import numpy as np

from clubb_jax.src.CLUBB_core.error_code import CLUBB_FATAL_ERROR
from clubb_jax.src.CLUBB_core import model_flags as _model_flags
from clubb_jax.src.CLUBB_core import ErrInfo
from clubb_jax.src.CLUBB_core.parameters_tunable import check_parameters

_iiPDF_ADG1 = _model_flags.iiPDF_ADG1
_iiPDF_ADG2 = _model_flags.iiPDF_ADG2
_iiPDF_3D_Luhar = _model_flags.iiPDF_3D_Luhar
_iiPDF_new = _model_flags.iiPDF_new
_iiPDF_TSDADG = _model_flags.iiPDF_TSDADG
_iiPDF_LY93 = _model_flags.iiPDF_LY93
_iiPDF_new_hybrid = _model_flags.iiPDF_new_hybrid
_saturation_bolton = _model_flags.saturation_bolton
_saturation_gfdl = _model_flags.saturation_gfdl
_saturation_flatau = _model_flags.saturation_flatau
_saturation_lookup = _model_flags.saturation_lookup
_ipdf_pre_advance_fields = _model_flags.ipdf_pre_advance_fields
_ipdf_pre_post_advance_fields = _model_flags.ipdf_pre_post_advance_fields
_l_explicit_turbulent_adv_wpxp = _model_flags.l_explicit_turbulent_adv_wpxp
_order_xm_wpxp = _model_flags.order_xm_wpxp
_order_xp2_xpyp = _model_flags.order_xp2_xpyp
_order_wp2_wp3 = _model_flags.order_wp2_wp3
_order_windm = _model_flags.order_windm


def parameterization_check(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    sclr_dim: int,
    edsclr_dim: int,
    thlm_forcing,
    rtm_forcing,
    um_forcing,
    vm_forcing,
    wm_zm,
    wm_zt,
    p_in_Pa,
    rho_zm,
    rho,
    exner,
    rho_ds_zm,
    rho_ds_zt,
    invrs_rho_ds_zm,
    invrs_rho_ds_zt,
    thv_ds_zm,
    thv_ds_zt,
    wpthlp_sfc,
    wprtp_sfc,
    upwp_sfc,
    vpwp_sfc,
    p_sfc,
    um,
    upwp,
    vm,
    vpwp,
    up2,
    vp2,
    rtm,
    wprtp,
    thlm,
    wpthlp,
    wp2,
    wp3,
    rtp2,
    thlp2,
    rtpthlp,
    prefix: str,
    wpsclrp_sfc,
    wpedsclrp_sfc,
    sclrm,
    wpsclrp,
    sclrp2,
    sclrprtp,
    sclrpthlp,
    sclrm_forcing,
    edsclrm,
    edsclrm_forcing,
    err_info: ErrInfo,
):
    """Determine what input variables may have NaN or invalid negative values.

    Description:
      This subroutine determines what input variables may have NaN values.
      In addition it checks to see if rho_zm, rho, exner, up2, vp2, rtm, thlm,
      wp2, rtp2, thlp2, or tau_zm have negative values.
    """
    del nzm, nzt, ngrdcol

    proc_name = "advance_clubb_core"
    operation = prefix + proc_name

    #-------- Input Nan Check ----------------------------------------------

    err_info = check_nan(thlm_forcing, "thlm_forcing", operation, err_info)
    err_info = check_nan(rtm_forcing, "rtm_forcing", operation, err_info)
    err_info = check_nan(um_forcing, "um_forcing", operation, err_info)
    err_info = check_nan(vm_forcing, "vm_forcing", operation, err_info)

    err_info = check_nan(wm_zm, "wm_zm", operation, err_info)
    err_info = check_nan(wm_zt, "wm_zt", operation, err_info)
    err_info = check_nan(p_in_Pa, "p_in_Pa", operation, err_info)
    err_info = check_nan(rho_zm, "rho_zm", operation, err_info)
    err_info = check_nan(rho, "rho", operation, err_info)
    err_info = check_nan(exner, "exner", operation, err_info)
    err_info = check_nan(rho_ds_zm, "rho_ds_zm", operation, err_info)
    err_info = check_nan(rho_ds_zt, "rho_ds_zt", operation, err_info)
    err_info = check_nan(invrs_rho_ds_zm, "invrs_rho_ds_zm", operation, err_info)
    err_info = check_nan(invrs_rho_ds_zt, "invrs_rho_ds_zt", operation, err_info)
    err_info = check_nan(thv_ds_zm, "thv_ds_zm", operation, err_info)
    err_info = check_nan(thv_ds_zt, "thv_ds_zt", operation, err_info)

    err_info = check_nan(um, "um", operation, err_info)
    err_info = check_nan(upwp, "upwp", operation, err_info)
    err_info = check_nan(vm, "vm", operation, err_info)
    err_info = check_nan(vpwp, "vpwp", operation, err_info)
    err_info = check_nan(up2, "up2", operation, err_info)
    err_info = check_nan(vp2, "vp2", operation, err_info)
    err_info = check_nan(rtm, "rtm", operation, err_info)
    err_info = check_nan(wprtp, "wprtp", operation, err_info)
    err_info = check_nan(thlm, "thlm", operation, err_info)
    err_info = check_nan(wpthlp, "wpthlp", operation, err_info)
    err_info = check_nan(wp2, "wp2", operation, err_info)
    err_info = check_nan(wp3, "wp3", operation, err_info)
    err_info = check_nan(rtp2, "rtp2", operation, err_info)
    err_info = check_nan(thlp2, "thlp2", operation, err_info)
    err_info = check_nan(rtpthlp, "rtpthlp", operation, err_info)

    err_info = check_nan(wpthlp_sfc, "wpthlp_sfc", operation, err_info)
    err_info = check_nan(wprtp_sfc, "wprtp_sfc", operation, err_info)
    err_info = check_nan(upwp_sfc, "upwp_sfc", operation, err_info)
    err_info = check_nan(vpwp_sfc, "vpwp_sfc", operation, err_info)
    err_info = check_nan(p_sfc, "p_sfc", operation, err_info)

    for sclr in range(sclr_dim):
        err_info = check_nan(
            sclrm_forcing[:, :, sclr], "sclrm_forcing", operation, err_info,
        )
        err_info = check_nan(
            wpsclrp_sfc[:, sclr], "wpsclrp_sfc", operation, err_info,
        )
        err_info = check_nan(sclrm[:, :, sclr], "sclrm", operation, err_info)
        err_info = check_nan(wpsclrp[:, :, sclr], "wpsclrp", operation, err_info)
        err_info = check_nan(sclrp2[:, :, sclr], "sclrp2", operation, err_info)
        err_info = check_nan(
            sclrprtp[:, :, sclr], "sclrprtp", operation, err_info,
        )
        err_info = check_nan(
            sclrpthlp[:, :, sclr], "sclrpthlp", operation, err_info,
        )

    for edsclr in range(edsclr_dim):
        err_info = check_nan(
            edsclrm_forcing[:, :, edsclr], "edsclrm_forcing", operation, err_info,
        )
        err_info = check_nan(
            wpedsclrp_sfc[:, edsclr], "wpedsclrp_sfc", operation, err_info,
        )
        err_info = check_nan(edsclrm[:, :, edsclr], "edsclrm", operation, err_info)

    #---------------------------------------------------------------------

    nan_fatal = err_info.any_fatal()

    err_info = check_negative(rtm, "rtm", operation, err_info)
    err_info = check_negative(p_in_Pa, "p_in_Pa", operation, err_info)
    err_info = check_negative(rho, "rho", operation, err_info)
    err_info = check_negative(rho_zm, "rho_zm", operation, err_info)
    err_info = check_negative(exner, "exner", operation, err_info)
    err_info = check_negative(rho_ds_zm, "rho_ds_zm", operation, err_info)
    err_info = check_negative(rho_ds_zt, "rho_ds_zt", operation, err_info)
    err_info = check_negative(
        invrs_rho_ds_zm, "invrs_rho_ds_zm", operation, err_info,
    )
    err_info = check_negative(
        invrs_rho_ds_zt, "invrs_rho_ds_zt", operation, err_info,
    )
    err_info = check_negative(thv_ds_zm, "thv_ds_zm", operation, err_info)
    err_info = check_negative(thv_ds_zt, "thv_ds_zt", operation, err_info)
    err_info = check_negative(up2, "up2", operation, err_info)
    err_info = check_negative(vp2, "vp2", operation, err_info)
    err_info = check_negative(wp2, "wp2", operation, err_info)
    err_info = check_negative(thlm, "thlm", operation, err_info)
    err_info = check_negative(rtp2, "rtp2", operation, err_info)
    err_info = check_negative(thlp2, "thlp2", operation, err_info)

    if prefix == "beginning of ":
        reset = err_info.reset_code()
        reset_mask = err_info.any_fatal() & jnp.logical_not(nan_fatal)
        err_info = err_info.replace(
            err_code=jnp.where(
                reset_mask,
                reset.err_code_or_default(),
                err_info.err_code_or_default(),
            ),
        )

    # Check the first levels for temperatures greater than 200K
    # Host diagnostic adaptation: retain the source's diagnostic-only test.
    def print_cold_levels(values):
        for i, k in np.argwhere(values < 190.0):
            print(
                "Liquid water potential temperature (thlm) < 190K ",
                "at grid column i = ", i + 1, " and grid level k = ", k + 1,
                ": thlm(", i + 1, ",", k + 1, ") = ", values[i, k],
                file=sys.stderr,
            )

    first_levels = thlm[:, :10]
    jax.lax.cond(
        jnp.any(first_levels < 190.0) & ~nan_fatal,
        lambda _: jax.debug.callback(print_cold_levels, first_levels, ordered=True),
        lambda _: None,
        operand=None,
    )

    return err_info


def check_negative(var, varname, operation, err_info: ErrInfo):
    """Checks for negative values in the var array and reports the column and
    vertical indices in which the negative values occur.

    """
    # Host I/O adaptation of Fortran's per-index writes, in level/column order.
    def print_negative(values, info):
        values = np.asarray(values)
        indices = np.argwhere(values.T < 0.0) if values.ndim == 2 else np.argwhere(values < 0.0)
        for index in indices:
            if isinstance(info, ErrInfo):
                print(info.global_error_header_host(), file=sys.stderr)
            k, i = ((int(index[0]) + 1, int(index[1]) + 1)
                    if values.ndim == 2 else (1, int(index[0]) + 1 if len(index) else 1))
            print(varname, " < 0 in ", operation,
                  " at i = ", i, " and k = ", k, file=sys.stderr)

    if not hasattr(err_info, "set_fatal"):
        if np.any(np.asarray(var) < 0.0):
            print_negative(np.asarray(var), None)
            err_info[...] = CLUBB_FATAL_ERROR
        return err_info
    jax.lax.cond(
        jnp.any(jnp.asarray(var) < 0.0),
        lambda _: jax.debug.callback(print_negative, var, err_info, ordered=True),
        lambda _: None,
        operand=None,
    )
    return err_info.set_fatal(mask=jnp.any(jnp.asarray(var) < 0.0))


def check_nan(var, varname, operation, err_info: ErrInfo):
    """Checks for a NaN in the var array and reports it.

    JAX treats any non-finite value as fatal, as does Fortran's IEEE check.
    """
    def print_nan(info):
        if isinstance(info, ErrInfo):
            print(info.global_error_header_host(), file=sys.stderr)
        print(varname, " is NaN in ", operation, file=sys.stderr)

    if not hasattr(err_info, "set_fatal"):
        if np.any(~np.isfinite(np.asarray(var))):
            print_nan(None)
            err_info[...] = CLUBB_FATAL_ERROR
        return err_info
    jax.lax.cond(
        jnp.any(~jnp.isfinite(jnp.asarray(var))),
        lambda _: jax.debug.callback(print_nan, err_info, ordered=True),
        lambda _: None,
        operand=None,
    )
    return err_info.set_fatal(
        mask=jnp.any(jnp.logical_not(jnp.isfinite(jnp.asarray(var))))
    )


def calculate_spurious_source(
    integral_after,
    integral_before,
    flux_top,
    flux_sfc,
    integral_forcing,
    dt,
):
    """Checks whether there is conservation within the column and returns any
    imbalance as spurious_source where spurious_source is defined negative
    for a spurious sink.
    """
    return (
        (jnp.asarray(integral_after) - jnp.asarray(integral_before))
        / jnp.asarray(dt)
        + jnp.asarray(flux_top)
        - jnp.asarray(flux_sfc)
        - jnp.asarray(integral_forcing)
    )


def _all_finite(*values):
    for value in values:
        if value is None:
            continue
        if not np.all(np.isfinite(np.asarray(value))):
            return False
    return True


def sfc_varnce_check(
    sclr_dim,
    wp2_sfc,
    up2_sfc,
    vp2_sfc,
    thlp2_sfc,
    rtp2_sfc,
    rtpthlp_sfc,
    sclrp2_sfc=None,
    sclrprtp_sfc=None,
    sclrpthlp_sfc=None,
):
    """Determine whether calc_surface_varnce outputs carry values that are nans."""
    values = [wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc]
    if int(sclr_dim) > 0:
        values.extend([sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc])
    for name, value in zip(
        ("wp2_sfc", "up2_sfc", "vp2_sfc", "thlp2_sfc", "rtp2_sfc", "rtpthlp_sfc",
         "sclrp2_sfc", "sclrprtp_sfc", "sclrpthlp_sfc"), values,
    ):
        if value is not None and np.any(~np.isfinite(np.asarray(value))):
            print(name, " is NaN in ", "calc_surface_varnce", file=sys.stderr)
    return _all_finite(*values)


def length_check(Lscale, Lscale_up, Lscale_down):
    """Report nonfinite combined/upward/downward lengths [m] and return host validity.

    The mixing-length JIT kernel checks finite values with JAX and invokes
    this diagnostic only inside its host error callback.
    """
    for name, value in (("Lscale", Lscale), ("Lscale_up", Lscale_up), ("Lscale_down", Lscale_down)):
        if value is not None and np.any(~np.isfinite(np.asarray(value))):
            print(name, " is NaN in ", "compute_mixing_length", file=sys.stderr)
    return _all_finite(Lscale, Lscale_up, Lscale_down)


def pdf_closure_check(closure_fields, pdf_params, sclr_dim=0, sclr_fields=None):
    """Determine whether pdf_closure output variables carry values that are NaNs."""
    values = []
    if isinstance(closure_fields, dict):
        values.extend(closure_fields.values())
    else:
        values.append(closure_fields)

    if hasattr(pdf_params, "_fields"):
        values.extend(getattr(pdf_params, field) for field in pdf_params._fields)
    else:
        values.append(pdf_params)

    if int(sclr_dim) > 0 and sclr_fields:
        values.extend(sclr_fields.values())
    return _all_finite(*values)


def rad_check(thlm, rcm, rtm, rim, cloud_frac, p_in_Pa, exner, rho_zm):
    """Checks radiation input variables. If they are < 0 it reports to the console."""
    values = [thlm, rcm, rtm, rim, cloud_frac, p_in_Pa, exner, rho_zm]
    if not _all_finite(*values):
        return False
    if any(np.any(np.asarray(value) < 0.0) for value in values):
        return False
    return bool(np.all(np.asarray(rtm) - np.asarray(rcm) >= 0.0))


def invalid_model_arrays(**arrays):
    """Checks for invalid floating point values in select model arrays.

    References:
      None
    """
    for name, value in arrays.items():
        if value is None:
            continue
        try:
            arr = np.asarray(value)
        except (TypeError, ValueError):
            continue
        if arr.dtype.kind not in "biufc":
            continue
        if not np.all(np.isfinite(arr)):
            print("NaN in", name, "model array", file=sys.stderr)
            return True
    return False


def check_clubb_settings(*, err_info, ngrdcol, edsclr_dim, sclr_idx, params, config_flags, **kwargs):
    """Subroutine to set up the model for execution.

    References:
      None

    Parameter validation is the first operation in current Fortran.  The
    l_do_expldiff_rtm_thlm eddy-scalar checks are ported; the remaining flag
    and mixed-setting checks are still handled by the JAX driver's
    supported-feature gates.
    """
    del kwargs
    # TODO: The remaining Fortran flag checks and their diagnostics are unported;
    # restoring that configuration machinery exceeds print-only work.
    err_info = check_parameters(ngrdcol, params, err_info)

    # The l_do_expldiff_rtm_thlm flag explicitly diffuses copies of thlm and
    # rtm that advance_windm_edsclrm stores in the eddy scalars
    # sclr_idx%iiedsclr_thl and sclr_idx%iiedsclr_rt.  With fewer than two
    # eddy scalars there is nowhere to store the copies, so the flag has no
    # effect.
    if config_flags.l_do_expldiff_rtm_thlm and edsclr_dim > 1:

        # The copies of thlm and rtm need distinct, valid eddy-scalar slots.
        if (sclr_idx.iiedsclr_thl < 1 or sclr_idx.iiedsclr_thl > edsclr_dim
                or sclr_idx.iiedsclr_rt < 1 or sclr_idx.iiedsclr_rt > edsclr_dim
                or sclr_idx.iiedsclr_thl == sclr_idx.iiedsclr_rt):
            print("l_do_expldiff_rtm_thlm = T requires distinct"
                  " iiedsclr_thl and iiedsclr_rt in 1:edsclr_dim.", file=sys.stderr)
            print("iiedsclr_thl = ", sclr_idx.iiedsclr_thl, file=sys.stderr)
            print("iiedsclr_rt = ", sclr_idx.iiedsclr_rt, file=sys.stderr)
            print("edsclr_dim = ", edsclr_dim, file=sys.stderr)
            print("Fatal error in check_clubb_settings_api", file=sys.stderr)
            # General error -> set all entries to clubb_fatal_error
            return err_info.set_fatal()

    return err_info


__all__ = [
    "CLUBB_FATAL_ERROR",
    "parameterization_check",
    "check_negative",
    "check_nan",
    "calculate_spurious_source",
    "sfc_varnce_check",
    "length_check",
    "pdf_closure_check",
    "rad_check",
    "invalid_model_arrays",
    "check_clubb_settings",
]
