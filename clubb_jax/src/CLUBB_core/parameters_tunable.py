"""Pure-Python implementation of selected ``parameters_tunable.F90`` APIs.

Description:
  This module contains tunable model parameters.  The purpose of the module is to make it
  easier for the clubb_tuner code to use the clubb_params vector without "knowing" any
  information about the individual parameters contained in the vector itself.
  It makes it easier to add
  new parameters to be tuned for, but does not make the CLUBB_core code itself any simpler.
  The parameters within the vector do not need to be the same variables used in the rest of
  CLUBB_core (see for e.g. nu1_vert_res_dep or lmin_coef).
  The parameters in the clubb_params vector only need to be those parameters for which we're not
  sure the correct value and we'd like to tune for.

References:
  None

Porting deviations:
The JAX driver needs ``get_param_names``, ``init_clubb_params``,
``calc_derrived_params``, and ``check_parameters``.  The Fortran tuner-only
``read_param_minmax`` and ``read_param_constraints`` routines are omitted.
Fortran ``pack_parameters`` and ``unpack_parameters`` are represented by
``PARAM_NAMES``/``PNAME_IDX`` plus a Python dictionary of defaults.  Fortran
mutates derived-type allocations in place; Python returns a ``NuVertResDep``
named tuple with NumPy arrays.
"""
from __future__ import annotations

import math
import sys

import numpy as np

from clubb_jax.src.Input_fields.namelist import read_namelist
from clubb_jax.src.CLUBB_core.nu_vert_res_dep import NuVertResDep


# ---------------------------------------------------------------------------
# These are referenced together often enough that it made sense to
# make a list of them.  Note that lmin_coef is the input parameter,
# while the actual lmin model constant is computed from this.
#***************************************************************
#                    ***** IMPORTANT *****
# If you change the order of the parameters in the parameter_indices,
# you will need to change the order of this list as well or the
# tuner will break!
#                    ***** IMPORTANT *****
#***************************************************************
# ---------------------------------------------------------------------------

PARAM_NAMES: list[str] = [
    "C1", "C1b", "C1c",
    "C2rt", "C2thl", "C2rtthl",
    "C4", "C_uu_shr", "C_uu_buoy",
    "C6rt", "C6rtb", "C6rtc",
    "C6thl", "C6thlb", "C6thlc",
    "C7", "C7b", "C7c",
    "C8", "C8b",
    "C10",
    "C11", "C11b", "C11c",
    "C12", "C13", "C14",
    "C_wp2_pr_dfsn", "C_wp3_pr_tp", "C_wp3_pr_turb", "C_wp3_pr_dfsn",
    "C_wp2_splat",
    "C6rt_Lscale0", "C6thl_Lscale0", "C7_Lscale0", "wpxp_L_thresh",
    "c_K", "c_K1", "nu1",
    "c_K2", "nu2",
    "c_K6", "nu6",
    "c_K8", "nu8",
    "c_K9", "nu9", "nu10",
    "c_K_hm", "c_K_hmb", "K_hm_min_coef", "nu_hm",
    "slope_coef_spread_DG_means_w", "pdf_component_stdev_factor_w",
    "coef_spread_DG_means_rt", "coef_spread_DG_means_thl",
    "gamma_coef", "gamma_coefb", "gamma_coefc",
    "mu", "beta", "lmin_coef",
    "omicron", "zeta_vrnce_rat", "upsilon_precip_frac_rat",
    "lambda0_stability_coef", "mult_coef",
    "taumin", "taumax",
    "Lscale_mu_coef", "Lscale_pert_coef",
    "alpha_corr", "Skw_denom_coef",
    "c_K10", "c_K10h",
    "thlp2_rad_coef", "thlp2_rad_cloud_frac_thresh",
    "up2_sfc_coef",
    "Skw_max_mag",
    "C_invrs_tau_bkgnd", "C_invrs_tau_sfc", "C_invrs_tau_shear",
    "C_invrs_tau_N2", "C_invrs_tau_N2_wp2", "C_invrs_tau_N2_xp2",
    "C_invrs_tau_N2_wpxp", "C_invrs_tau_N2_clear_wp3",
    "C_invrs_tau_wpxp_Ri", "C_invrs_tau_wpxp_N2_thresh",
    "xp3_coef_base", "xp3_coef_slope",
    "altitude_threshold", "rtp2_clip_coef",
    "Cx_min", "Cx_max",
    "Richardson_num_min", "Richardson_num_max",
    "a3_coef_min", "a_const", "bv_efold", "wpxp_Ri_exp", "z_displace",
]
assert len(PARAM_NAMES) == 102, f"Expected 102 params, got {len(PARAM_NAMES)}"

# Case-insensitive lookup: lower-case name → 0-based index
_NAME_TO_IDX: dict[str, int] = {n.lower(): i for i, n in enumerate(PARAM_NAMES)}

# Case-sensitive name → 0-based index (the single source of truth for parameter ordering; numerical_check
# imports this for check_clubb_settings).
PNAME_IDX: dict[str, int] = {n: i for i, n in enumerate(PARAM_NAMES)}

# Sentinel for a fatal validation error (error_code.F90); used by check_parameters below.
CLUBB_FATAL_ERROR = 99
_EPS64 = np.finfo(np.float64).eps
_TINY64 = np.finfo(np.float64).tiny
_HUGE64 = np.finfo(np.float64).max

# Hard rejection bounds from parameters_tunable.F90.  Most parameters are
# nonnegative with no finite upper bound; only the source exceptions are
# overridden below.  Rows are [minimum, maximum].
PARAMETER_HARD_BOUNDS = np.vstack(
    (
        np.zeros(len(PARAM_NAMES), dtype=np.float64),
        np.full(len(PARAM_NAMES), _HUGE64, dtype=np.float64),
    )
)
for _name in ("C_uu_shr", "C_uu_buoy", "C7", "C7b", "C11", "C11b",
              "upsilon_precip_frac_rat"):
    PARAMETER_HARD_BOUNDS[1, PNAME_IDX[_name]] = 1.0
for _name in ("slope_coef_spread_DG_means_w", "pdf_component_stdev_factor_w"):
    PARAMETER_HARD_BOUNDS[0, PNAME_IDX[_name]] = _TINY64
for _name in ("coef_spread_DG_means_rt", "coef_spread_DG_means_thl"):
    PARAMETER_HARD_BOUNDS[1, PNAME_IDX[_name]] = 1.0 - _EPS64
PARAMETER_HARD_BOUNDS[1, PNAME_IDX["beta"]] = 3.0
PARAMETER_HARD_BOUNDS[:, PNAME_IDX["omicron"]] = (_TINY64, 1.0)
PARAMETER_HARD_BOUNDS[0, PNAME_IDX["zeta_vrnce_rat"]] = -1.0 + _EPS64
PARAMETER_HARD_BOUNDS[:, PNAME_IDX["a3_coef_min"]] = (1.0, 3.0)

# ---------------------------------------------------------------------------
# Default parameter values (mirrors set_default_parameters in Fortran)
# ---------------------------------------------------------------------------

_DEFAULTS: dict[str, float] = {
    "C1": 1.0, "C1b": 1.0, "C1c": 1.0,
    "C2rt": 2.0, "C2thl": 2.0, "C2rtthl": 2.0,
    "C4": 2.0, "C_uu_shr": 0.4, "C_uu_buoy": 0.3,
    "C6rt": 2.0, "C6rtb": 2.0, "C6rtc": 1.0,
    "C6thl": 2.0, "C6thlb": 2.0, "C6thlc": 1.0,
    "C7": 0.5, "C7b": 0.5, "C7c": 0.5,
    "C8": 0.5, "C8b": 0.02,
    "C10": 3.3,
    "C11": 0.4, "C11b": 0.4, "C11c": 0.5,
    "C12": 1.0, "C13": 0.1, "C14": 1.0,
    "C_wp2_pr_dfsn": 0.0, "C_wp3_pr_tp": 0.0, "C_wp3_pr_turb": 0.0,
    "C_wp3_pr_dfsn": 0.0, "C_wp2_splat": 2.0,
    "C6rt_Lscale0": 14.0, "C6thl_Lscale0": 14.0, "C7_Lscale0": 0.85,
    "wpxp_L_thresh": 60.0,
    "c_K": 0.2, "c_K1": 0.2,
    "nu1": 20.0, "c_K2": 0.025, "nu2": 1.0,
    "c_K6": 0.375, "nu6": 5.0,
    "c_K8": 5.0, "nu8": 20.0,
    "c_K9": 0.1, "nu9": 10.0, "nu10": 0.0,
    "c_K_hm": 0.75, "c_K_hmb": 0.75, "K_hm_min_coef": 0.1, "nu_hm": 1.5,
    "slope_coef_spread_DG_means_w": 21.0, "pdf_component_stdev_factor_w": 1.0,
    "coef_spread_DG_means_rt": 0.8, "coef_spread_DG_means_thl": 0.8,
    "gamma_coef": 0.25, "gamma_coefb": 0.25, "gamma_coefc": 5.0,
    "mu": 1.0e-3, "beta": 1.0, "lmin_coef": 0.5,
    "omicron": 0.5, "zeta_vrnce_rat": 0.0, "upsilon_precip_frac_rat": 0.55,
    "lambda0_stability_coef": 0.03, "mult_coef": 0.5,
    "taumin": 90.0, "taumax": 3600.0,
    "Lscale_mu_coef": 2.0, "Lscale_pert_coef": 0.1,
    "alpha_corr": 0.15, "Skw_denom_coef": 4.0,
    "c_K10": 1.0, "c_K10h": 1.0,
    "thlp2_rad_coef": 1.0, "thlp2_rad_cloud_frac_thresh": 0.1,
    "up2_sfc_coef": 4.0,
    "Skw_max_mag": 10.0,
    "C_invrs_tau_bkgnd": 1.1, "C_invrs_tau_sfc": 0.1, "C_invrs_tau_shear": 0.15,
    "C_invrs_tau_N2": 0.4, "C_invrs_tau_N2_wp2": 0.2, "C_invrs_tau_N2_xp2": 0.05,
    "C_invrs_tau_N2_wpxp": 0.0, "C_invrs_tau_N2_clear_wp3": 1.0,
    "C_invrs_tau_wpxp_Ri": 0.35, "C_invrs_tau_wpxp_N2_thresh": 3.3e-4,
    "xp3_coef_base": 0.25, "xp3_coef_slope": 0.01,
    "altitude_threshold": 100.0, "rtp2_clip_coef": 0.5,
    "Cx_min": 0.33, "Cx_max": 0.95,
    "Richardson_num_min": 0.25, "Richardson_num_max": 400.0,
    "a3_coef_min": 1.0, "a_const": 1.8, "bv_efold": 5.0,
    "wpxp_Ri_exp": 0.5, "z_displace": 25.0,
}
assert set(_DEFAULTS.keys()) == set(PARAM_NAMES), "Defaults must cover all 102 parameters"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_param_names() -> list[str]:
    """Return the 102 parameter names in canonical Fortran ordering."""
    return list(PARAM_NAMES)


def init_clubb_params(ngrdcol: int, filename: str) -> np.ndarray:
    """Read a namelist containing the model parameters.

    Description:
    Read a namelist containing the model parameters

    References:
    None
    """
    nparams = 102

    # Set the default tunable parameter values
    values = np.tile(
        np.array([_DEFAULTS[n] for n in PARAM_NAMES], dtype=np.float64),
        (ngrdcol, 1),
    )

    # If the filename is empty, assume we're using a `working' set of
    # parameters that are set statically here (handy for host models).
    # Read the namelist
    nml = read_namelist(filename)
    for raw_name, val in nml.items():
        key = raw_name.lower().strip()
        if key in _NAME_TO_IDX:
            idx = _NAME_TO_IDX[key]
            arr = np.asarray(val, dtype=np.float64)
            if arr.ndim == 0:
                values[:, idx] = float(arr)
            else:
                arr = arr.ravel()
                if arr.size == 1:
                    values[:, idx] = float(arr[0])
                elif arr.size == ngrdcol:
                    values[:, idx] = arr
                else:
                    raise ValueError(
                        f"{raw_name} must be scalar or have ngrdcol={ngrdcol} values; "
                        f"got {arr.size}."
                    )

    return values


def calc_derrived_params(
    gr,
    ngrdcol: int,
    grid_type: int,
    deltaz: np.ndarray,
    clubb_params: np.ndarray,
    l_prescribed_avg_deltaz: bool,
) -> tuple[NuVertResDep, float, float]:
    """Calculate parameters that should be derrived from other quantities.

    Description:
      Calculates clubb parameters that should be derrived from other quantities.

      Adjusts the values of background eddy diffusivity based on
      vertical grid spacing.
      This code was made into a public subroutine so that it may be
      called multiple times per model run in scenarios where grid
      altitudes, and hence average grid spacing, change through space
      and/or time.  This occurs, for example, when CLUBB is
      implemented in WRF.  --ldgrant Jul 2010
    """
    #------------------------------ Constant Parameters ------------------------------
    # Fixed value for minimum value for the length scale.
    lmin_deltaz = 40.0

    # It was decided after some experimentation, that the best
    # way to produce grid independent results is to set lmin to be
    # some fixed value. -dschanen 21 May 2007
    # TODO: using "clubb_params(ngrdcol,ilmin_coef)", but lmin should really be
    # changed to dimension(ngrdcol) to avoid this
    lmin = float(clubb_params[ngrdcol - 1, 61]) * lmin_deltaz
    if lmin < 1.0:
        raise ValueError("lmin is < 1.0")

    # Using ngrdcol here as well for temporary backward compatibility, same as above
    Skw_max = float(clubb_params[ngrdcol - 1, 78])
    inner = 4.0 * (1.0 - 0.4) ** 3 + Skw_max ** 2
    mixt_frac_max_mag = 1.0 - 0.5 * (1.0 - Skw_max / math.sqrt(inner))
    # Known magic number

    #------------------------------ Local Variables ------------------------------
    # Average grid box height   [m]
    deltaz = np.asarray(deltaz, dtype=np.float64).ravel()
    avg_deltaz = np.empty(ngrdcol, dtype=np.float64)

    if l_prescribed_avg_deltaz or grid_type == 1:
        avg_deltaz[:] = deltaz
    elif grid_type == 2:
        # Stretched (unevenly-spaced) grid:  stretched thermodynamic level
        # input.
        # Find the average deltaz over the stretched grid based on
        # thermodynamic level inputs.
        zt = np.asarray(gr.zt)  # (ngrdcol, nzt)
        for i in range(ngrdcol):
            avg_deltaz[i] = (zt[i, -1] - zt[i, 0]) / max(1, zt.shape[1] - 1)
    elif grid_type == 3:
        # CLUBB is implemented in a host model, or is using grid_type = 3
        # Find the average deltaz over the grid based on momentum level
        # inputs.
        zm = np.asarray(gr.zm)  # (ngrdcol, nzm)
        for i in range(ngrdcol):
            avg_deltaz[i] = (zm[i, -1] - zm[i, 0]) / max(1, zm.shape[1] - 1)
    else:
        avg_deltaz[:] = deltaz

    # Flag for adjusting the values of the constant background eddy diffusivity
    # coefficients based on the average vertical grid spacing.  If this flag is
    # turned off, the values of the various nu coefficients will remain as they
    # are declared in the tunable_parameters.in file.

    # The size of the average vertical grid spacing that serves as a threshold
    # for when to increase the size of the background eddy diffusivity
    # coefficients (nus) by a certain factor above what the background
    # coefficients are specified to be in tunable_parameters.in.  At any average
    # grid spacing at or below this value, the values of the background
    # diffusivities remain the same.  However, at any average vertical grid
    # spacing above this value, the values of the background eddy diffusivities
    # are increased.  Traditionally, the threshold grid spacing has been set to
    # 40.0 meters.  This is only relevant if l_adj_low_res_nu is turned on.
    grid_spacing_thresh = 40.0

    # The factor by which to multiply the coefficients of background eddy
    # diffusivity if the grid spacing threshold is exceeded and l_adj_low_res_nu
    # is turned on.
    mult_factor = np.ones(ngrdcol, dtype=np.float64)
    for i in range(ngrdcol):
        if avg_deltaz[i] > grid_spacing_thresh:
            mult_factor[i] = 1.0 + float(clubb_params[i, 66]) * math.log(
                avg_deltaz[i] / grid_spacing_thresh
            )

    # The nu's are chosen for deltaz <= 40 m. Looks like they must
    # be adjusted for larger grid spacings (Vince Larson)

    # Use a constant mult_factor so nu does not depend on grid spacing
    nu1   = clubb_params[:, 38] * mult_factor
    nu2   = clubb_params[:, 40] * mult_factor
    nu6   = clubb_params[:, 42] * mult_factor
    nu8   = clubb_params[:, 44] * mult_factor   # zt-level
    nu9   = clubb_params[:, 46] * mult_factor
    nu10  = clubb_params[:, 47] * mult_factor   # zt-level (disabled in ARM)
    nu_hm = clubb_params[:, 51] * mult_factor   # zt-level

    nu_vert_res_dep = NuVertResDep(
        nzm=int(gr.nzm),
        nu1=nu1.copy(),
        nu2=nu2.copy(),
        nu6=nu6.copy(),
        nu8=nu8.copy(),
        nu9=nu9.copy(),
        nu10=nu10.copy(),
        nu_hm=nu_hm.copy(),
    )

    return nu_vert_res_dep, lmin, mixt_frac_max_mag


def get_parameter_hard_bounds(nparams_in: int = len(PARAM_NAMES)) -> np.ndarray:
    """Return the authoritative 2-by-nparams hard configuration bounds."""
    if nparams_in != len(PARAM_NAMES):
        raise ValueError(
            f"get_parameter_hard_bounds: expected {len(PARAM_NAMES)} parameters, "
            f"got {nparams_in}"
        )
    return PARAMETER_HARD_BOUNDS.copy()


def check_parameter_hard_bounds(
    ngrdcol: int,
    clubb_params: np.ndarray,
    err_info,
):
    """Apply the same inclusive hard-bound table as Fortran."""
    values = np.asarray(clubb_params, dtype=np.float64)
    if values.shape != (ngrdcol, len(PARAM_NAMES)):
        raise ValueError(
            f"clubb_params must have shape ({ngrdcol}, {len(PARAM_NAMES)}); "
            f"got {values.shape}"
        )

    existing_codes = (
        err_info.err_code_or_default()
        if hasattr(err_info, "err_code_or_default")
        else err_info.err_code
    )
    err_code = np.asarray(existing_codes, dtype=np.int32).copy()
    below = values < PARAMETER_HARD_BOUNDS[0][None, :]
    above = values > PARAMETER_HARD_BOUNDS[1][None, :]
    invalid = below | above
    for i, k in np.argwhere(invalid):
        name = PARAM_NAMES[int(k)]
        lower = PARAMETER_HARD_BOUNDS[0, k]
        upper = PARAMETER_HARD_BOUNDS[1, k]
        print(
            f"{name} = {values[i, k]} is outside hard configuration bounds: "
            f"{lower} {upper}",
            file=sys.stderr,
        )
        err_code[int(i)] = CLUBB_FATAL_ERROR
    return err_info._replace(err_code=err_code)


def check_parameters(
    ngrdcol: int,
    clubb_params: np.ndarray,
    err_info,
):
    """Validate Fortran's hard bounds and required paired equalities."""
    values = np.asarray(clubb_params, dtype=np.float64)
    err_info = check_parameter_hard_bounds(ngrdcol, values, err_info)
    err_code = np.asarray(err_info.err_code, dtype=np.int32).copy()

    equal_pairs = (
        ("C6rt", "C6thl"),
        ("C6rtb", "C6thlb"),
        ("C6rtc", "C6thlc"),
        ("C6rt_Lscale0", "C6thl_Lscale0"),
    )
    for i in range(ngrdcol):
        for left_name, right_name in equal_pairs:
            left = values[i, PNAME_IDX[left_name]]
            right = values[i, PNAME_IDX[right_name]]
            tolerance = abs(left + right) * 0.5 * _EPS64
            if abs(left - right) > tolerance:
                print(
                    f"{left_name} = {left}, {right_name} = {right} "
                    f"({left_name} and {right_name} must be equal)",
                    file=sys.stderr,
                )
                err_code[i] = CLUBB_FATAL_ERROR

    return err_info._replace(err_code=err_code)
