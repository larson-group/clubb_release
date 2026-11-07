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
    """This subroutine determines what input variables may have NaN values.
    In addition it checks to see if rho_zm, rho, exner, up2, vp2, rtm, thlm,
    wp2, rtp2, thlp2, or tau_zm have negative values.

    Returns:
        err_info

    Fortran declarations (parameterization_check in src/CLUBB_core/numerical_check.F90):

        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol, &
          sclr_dim, &
          edsclr_dim

        ! Input variables
        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
          thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
          rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
          um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
          vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
          wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
          p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
          rho,             & ! Air density on thermodynamic levels       [kg/m^3]
          exner,           & ! Exner function (thermodynamic levels)     [-]
          rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
          invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
          thv_ds_zt        ! Dry, base-state theta_v on thermo. levs.  [K]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          wm_zm,           & ! w mean wind component on momentum levels  [m/s]
          rho_zm,          & ! Air density on momentum levels            [kg/m^3]
          rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
          invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
          thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
          wpthlp_sfc,   & ! w' theta_l' at surface.   [(m K)/s]
          wprtp_sfc,    & ! w' r_t' at surface.       [(kg m)/( kg s)]
          upwp_sfc,     & ! u'w' at surface.          [m^2/s^2]
          vpwp_sfc,     & ! v'w' at surface.          [m^2/s^2]
          p_sfc           ! pressure at surface       [Pa]

        ! These are prognostic or are planned to be in the future
        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
          um,      & ! u mean wind component (thermodynamic levels)   [m/s]
          vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
          rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
          thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
          wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
          upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
          vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
          up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
          vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
          wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
          wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
          rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
          thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
          rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
          wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

        character(len=*), intent(in) :: prefix ! Location where subroutine is called

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,sclr_dim) :: &
          wpsclrp_sfc    ! Scalar flux at surface [units m/s]

        real( kind = core_rknd ), intent(in), dimension(ngrdcol,edsclr_dim) :: &
          wpedsclrp_sfc ! Eddy-Scalar flux at surface      [units m/s]

        real( kind = core_rknd ), intent(in),dimension(ngrdcol,nzt,sclr_dim) :: &
          sclrm,         & ! Passive scalar mean      [units vary]
          sclrm_forcing    ! Passive scalar forcing   [units / s]

        real( kind = core_rknd ), intent(in),dimension(ngrdcol,nzm,sclr_dim) :: &
          wpsclrp,       & ! w'sclr'                  [units vary]
          sclrp2,        & ! sclr'^2                  [units vary]
          sclrprtp,      & ! sclr'rt'                 [units vary]
          sclrpthlp        ! sclr'thl'                [units vary]

        real( kind = core_rknd ), intent(in),dimension(ngrdcol,nzt,edsclr_dim) :: &
          edsclrm,         & ! Eddy passive scalar mean    [units vary]
          edsclrm_forcing    ! Eddy passive scalar forcing [units / s]

        ! Input/Output variables
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
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
    """Checks for negative values in the var array and reports
    the column and vertical indices in which the negative values occur.

    Fortran declarations (check_negative in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: var(:,:)

        character(len=*), intent(in):: &
          varname,   & ! Variable being examined
          operation    ! Procedure calling check_zero

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
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
        jnp.any(var < 0.0),
        lambda _: jax.debug.callback(print_negative, var, err_info, ordered=True),
        lambda _: None,
        operand=None,
    )
    return err_info.set_fatal(mask=jnp.any(var < 0.0))


def check_nan(var, varname, operation, err_info: ErrInfo):
    """Checks for a NaN in the var array and reports it.

    Fortran declarations (check_nan_2d in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in), dimension(:,:) :: var ! Variable being examined

        character(len=*), intent(in):: &
          varname,  & ! Name of variable
          operation   ! Procedure calling check_nan

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
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
        jnp.any(~jnp.isfinite(var)),
        lambda _: jax.debug.callback(print_nan, err_info, ordered=True),
        lambda _: None,
        operand=None,
    )
    return err_info.set_fatal(
        mask=jnp.any(jnp.logical_not(jnp.isfinite(var)))
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

    Returns:
        spurious_source

    Fortran declarations (calculate_spurious_source in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          integral_after, &   ! Vertically-integrated quantity after dt time  [units vary]
          integral_before, &  ! Vertically-integrated quantity before dt time [units vary]
          flux_top, &         ! Total flux at the top of the domain           [units vary]
          flux_sfc, &         ! Total flux at the bottom of the domain        [units vary]
          integral_forcing, & ! Vertically-integrated forcing                 [units vary]
          dt                  ! Timestep size                                 [s]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: spurious_source ! [units vary]
    """
    spurious_source = (
        (integral_after - integral_before)
        / dt
        + flux_top
        - flux_sfc
        - integral_forcing
    )
    return spurious_source


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
    """This subroutine determines if any of the output
           variables for the calc_surface_varnce subroutine carry values that
           are nans.

           Joshua Fasching February 2008

    Fortran declarations (sfc_varnce_check in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          sclr_dim

        real( kind = core_rknd ),intent(in) :: &
          wp2_sfc,     & ! Vertical velocity variance        [m^2/s^2]
          up2_sfc,     & ! u'^2                              [m^2/s^2]
          vp2_sfc,     & ! u'^2                              [m^2/s^2]
          thlp2_sfc,   & ! thetal variance                   [K^2]
          rtp2_sfc,    & ! rt variance                       [(kg/kg)^2]
          rtpthlp_sfc    ! thetal rt covariance              [kg K/kg]

        real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
          sclrp2_sfc,    & ! Passive scalar variance                 [units^2]
          sclrprtp_sfc,  & ! Passive scalar r_t covariance           [units kg/kg]
          sclrpthlp_sfc    ! Passive scalar theta_l covariance       [units K]

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
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
    """This subroutine determines if any of the output
            variables for the length_new subroutine carry values that
            are NaNs.

            Joshua Fasching February 2008

    The JIT-compiled mixing-length kernel checks for finite values with JAX and
    calls this diagnostic only inside its host error callback.

    Fortran declarations (length_check in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzt

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), dimension(nzt), intent(in) :: &
          Lscale,     & ! Mixing length                 [m]
          Lscale_up,  & ! Upward mixing length          [m]
          Lscale_down   ! Downward mixing length        [m]

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    for name, value in (("Lscale", Lscale), ("Lscale_up", Lscale_up), ("Lscale_down", Lscale_down)):
        if value is not None and np.any(~np.isfinite(np.asarray(value))):
            print(name, " is NaN in ", "compute_mixing_length", file=sys.stderr)
    return _all_finite(Lscale, Lscale_up, Lscale_down)


def pdf_closure_check(closure_fields, pdf_params, sclr_dim=0, sclr_fields=None):
    """This subroutine determines if any of the output
       variables for the pdf_closure subroutine carry values that
       are NaNs.

     Joshua Fasching February 2008

    Fortran declarations (pdf_closure_check in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz, &
          sclr_dim

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), dimension(nz), intent(in) :: &
          wp4,             & ! w'^4                  [m^4/s^4]
          wprtp2,          & ! w' r_t'               [(m kg)/(s kg)]
          wp2rtp,          & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
          wpthlp2,         & ! w' th_l'^2            [(m K^2)/s]
          wp2thlp,         & ! w'^2 th_l'            [(m^2 K)/s^2]
          cloud_frac,      & ! Cloud fraction        [-]
          rcm,             & ! Mean liquid water     [kg/kg]
          wpthvp,          & ! Buoyancy flux         [(K m)/s]
          wp2thvp,         & ! w'^2 th_v'            [(m^2 K)/s^2]
          wp2up,           & ! w'^2 u'               [m^3/s^3]
          rtpthvp,         & ! r_t' th_v'            [(kg K)/kg]
          thlpthvp,        & ! th_l' th_v'           [K^2]
          wprcp,           & ! w' r_c'               [(m kg)/(s kg)]
          wp2rcp,          & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
          rtprcp,          & ! r_t' r_c'             [(kg^2)/(kg^2)]
          thlprcp,         & ! th_l' r_c'            [(K kg)/kg]
          rcp2,            & ! r_c'^2                [(kg^2)/(kg^2)]
          wprtpthlp,       & ! w' r_t' th_l'         [(m kg K)/(s kg)]
          crt_1, crt_2,  &
          cthl_1, cthl_2

        type(pdf_parameter), intent(in) :: &
          pdf_params        ! PDF parameters          [units vary]

        ! Input (Optional passive scalar variables)
        real( kind = core_rknd ), dimension(nz,sclr_dim), intent(in) :: &
          sclrpthvp, &
          sclrprcp, &
          wpsclrp2, &
          wpsclrprtp, &
          wpsclrpthlp, &
          wp2sclrp

        type(stats_type), intent(in) :: &
          stats

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
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
    """Checks radiation input variables. If they are < 0 it reports
    to the console.

    Fortran declarations (rad_check in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, &
          nzm, &
          nzt

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          thlm,           & ! Liquid Water Potential Temperature   [K/s]
          rcm,            & ! Liquid Water Mixing Ratio            [kg/kg]
          rtm,            & ! Total Water Mixing Ratio             [kg/kg]
          rim,            & ! Ice Water Mixing Ratio               [kg/kg]
          cloud_frac,     & ! Cloud Fraction                       [-]
          p_in_Pa,        & ! Pressure                             [Pa]
          exner             ! Exner Function                       [-]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rho_zm            ! Air Density                          [kg/m^3]

        !--------------------------- Input/Output Variables ---------------------------
        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    values = [thlm, rcm, rtm, rim, cloud_frac, p_in_Pa, exner, rho_zm]
    if not _all_finite(*values):
        return False
    if any(np.any(np.asarray(value) < 0.0) for value in values):
        return False
    return bool(np.all(np.asarray(rtm) - np.asarray(rcm) >= 0.0))


def invalid_model_arrays(**arrays):
    """Checks for invalid floating point values in select model arrays.

    Returns:
        invalid_model_arrays

    Fortran declarations (invalid_model_arrays in src/CLUBB_core/numerical_check.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          hydromet_dim, &
          sclr_dim, &
          edsclr_dim

        character(len=10), dimension(hydromet_dim), intent(in) :: &
          hydromet_list

        real( kind = core_rknd ), dimension(nzt), intent(in) :: &
          um,       & ! eastward grid-mean wind comp. (thermo. levs.)  [m/s]
          vm,       & ! northward grid-mean wind comp. (thermo. levs.) [m/s]
          rtm,      & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
          thlm,     & ! liq. water pot. temp., th_l (thermo. levels)   [K]
          wp3,      & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
          wp2thvp,  & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
          wp2up       ! < w'^2 u' > (thermodynamic levels)             [m^3/s^3]

        real( kind = core_rknd ), dimension(nzm), intent(in) :: &
          wprtp,    & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
          wpthlp,   & ! w'th_l' (momentum levels)                      [(m/s) K]
          rtp2,     & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
          thlp2,    & ! th_l'^2 (momentum levels)                      [K^2]
          rtpthlp,  & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
          wp2,      & ! w'^2 (momentum levels)                         [m^2/s^2]
          rtpthvp,  & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
          thlpthvp    ! < th_l' th_v' > (momentum levels)              [K^2]

        real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
          hydromet    ! Array of hydrometeors                          [units vary]

        real( kind = core_rknd ), dimension(nzt,sclr_dim), intent(in) :: &
          sclrm    ! Passive scalar mean (thermo. levels)              [units vary]

        real( kind = core_rknd ), dimension(nzt,edsclr_dim), intent(in) :: &
          edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]
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
