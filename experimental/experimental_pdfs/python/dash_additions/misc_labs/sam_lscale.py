"""Native CLUBB direct-Lscale diagnostic driven by resolved raw-SAM moments.

The interactive Misc laboratory must not duplicate the parcel/mixing-length algorithm in
Python.  This module only prepares the complete-column inputs that are
available from SAM micro files, then calls the existing F2PY wrapper around
``mixing_length.calc_Lscale_directly``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import sys
from threading import Lock

import numpy as np

from utilities.sam_3d_reference import ProfileMoments, load_profile_moments


_REPO_ROOT = Path(__file__).resolve().parents[2]
_PYTHON_API_ROOT = _REPO_ROOT / "clubb_python_api"
_FORTRAN_API_LOCK = Lock()


class SamLscaleError(RuntimeError):
    """Raised when the optional native CLUBB Lscale diagnostic is unavailable."""


@dataclass(frozen=True)
class SamLscaleProfile:
    """Direct-Lscale result and explicit provenance/approximations."""

    source_file: str
    elapsed_seconds: int
    elapsed_minutes: float
    height_m: np.ndarray
    lscale_m: np.ndarray
    lscale_up_m: np.ndarray
    lscale_down_m: np.ndarray
    tke_proxy_m2_s2: np.ndarray
    lmin_m: float
    newmu_m_inv: float
    path_kind: str
    approximation: str


def _api():
    """Import the compiled API lazily so an unavailable binding cannot break Dash."""

    if str(_PYTHON_API_ROOT) not in sys.path:
        sys.path.insert(0, str(_PYTHON_API_ROOT))
    # ``clubb_python`` is source-controlled, while its F2PY extension belongs
    # to one selected CMake build tree.  Dash is intentionally not installed
    # into that build tree, so discover the existing runtime artifact instead
    # of making each user hand-edit PYTHONPATH before opening Misc.
    runtime_roots = sorted(
        _REPO_ROOT.glob("build/*_PYTHON/clubb_python_api/f2py_runtime"),
        key=lambda path: ("intel" not in path.as_posix().lower(), path.as_posix()),
    )
    for runtime_root in reversed(runtime_roots):
        if str(runtime_root) not in sys.path:
            sys.path.insert(0, str(runtime_root))
    try:
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except (ImportError, OSError) as error:
        raise SamLscaleError(
            "The compiled CLUBB Python API is unavailable; build it with "
            "./compile.py -python before using the SAM Lscale diagnostic."
        ) from error
    return clubb_api, ErrInfo


def _momentum_heights(thermodynamic_heights: np.ndarray) -> np.ndarray:
    """Infer the one-more-level momentum grid around SAM's cell-center heights."""

    zt = np.asarray(thermodynamic_heights, dtype=float)
    if zt.size < 2 or not np.all(np.diff(zt) > 0.0):
        raise SamLscaleError("SAM thermodynamic heights must be strictly increasing.")
    zm = np.empty(zt.size + 1, dtype=float)
    zm[1:-1] = 0.5 * (zt[:-1] + zt[1:])
    zm[0] = zt[0] - 0.5 * (zt[1] - zt[0])
    zm[-1] = zt[-1] + 0.5 * (zt[-1] - zt[-2])
    return zm


def _exner(pressure_pa: np.ndarray) -> np.ndarray:
    """Use CLUBB's standard dry-air Exner definition from raw SAM pressure."""

    pressure = np.maximum(np.asarray(pressure_pa, dtype=float), 1.0)
    return (pressure / 100000.0) ** (287.04 / 1004.67)


def _parameter_value(api, params: np.ndarray, name: str) -> float:
    names = [str(item).strip().lower() for item in api.get_param_names()]
    try:
        return float(params[0, names.index(name.lower())])
    except ValueError as error:
        raise SamLscaleError(f"The compiled Python API does not expose CLUBB parameter {name!r}.") from error


def _prepare_inputs(profile: ProfileMoments) -> dict[str, np.ndarray]:
    """Map resolved SAM moments to the native direct-Lscale inputs.

    SAM micro output lacks horizontal velocity.  The only turbulence-energy
    estimate possible from this dataset is the documented isotropic proxy
    ``em = 3/2 w'^2``.  Every thermodynamic field and covariance is otherwise
    a direct resolved horizontal reduction.
    """

    zt = np.asarray(profile.height_m, dtype=float)
    zm = _momentum_heights(zt)
    covariance = np.asarray(profile.covariance, dtype=float)
    return {
        "zt": zt,
        "zm": zm,
        "rtm": np.asarray(profile.mean[:, 1], dtype=float),
        "thlm": np.asarray(profile.mean[:, 2], dtype=float),
        "rcm": np.asarray(profile.mean_rc, dtype=float),
        "rtp2": np.maximum(covariance[:, 1, 1], 1.0e-18),
        "thlp2": np.maximum(covariance[:, 2, 2], 1.0e-12),
        "rtpthlp": covariance[:, 1, 2],
        "em": 1.5
        * np.maximum(
            np.interp(zm, zt, np.maximum(covariance[:, 0, 0], 0.0)), 0.0
        ),
        "pressure": np.asarray(profile.pressure_pa, dtype=float),
    }


def _compute_native_lscale_profile(
    profile: ProfileMoments,
    *,
    rt_path: np.ndarray,
    thl_path: np.ndarray,
    path_kind: str,
    approximation: str,
) -> SamLscaleProfile:
    """Call CLUBB's direct parcel calculation for one supplied column path.

    The lock is essential: the F2PY wrapper stores CLUBB derived types in a
    process-global Fortran module.  One Misc callback therefore cannot change
    the native grid/PDF state while another is inside the same diagnostic.
    """

    inputs = _prepare_inputs(profile)
    rt_path = np.asarray(rt_path, dtype=float)
    thl_path = np.asarray(thl_path, dtype=float)
    expected_shape = inputs["rtm"].shape
    if rt_path.shape != expected_shape or thl_path.shape != expected_shape:
        raise SamLscaleError(
            "The supplied thermodynamic path must have one r_t and theta_l value "
            f"per SAM level ({expected_shape[0]} values)."
        )
    if not np.all(np.isfinite(rt_path)) or not np.all(np.isfinite(thl_path)):
        raise SamLscaleError("The supplied thermodynamic path contains non-finite values.")

    api, ErrInfo = _api()
    with _FORTRAN_API_LOCK:
        ngrdcol = 1
        nzt = int(inputs["zt"].size)
        nzm = int(inputs["zm"].size)
        err_info = ErrInfo(ngrdcol=ngrdcol)
        api.init_err_info(ngrdcol)
        flags = api.get_default_config_flags()
        api.init_config_flags(flags)
        grid, _ = api.setup_grid(
            nzmax=nzm,
            ngrdcol=ngrdcol,
            sfc_elevation=np.zeros(ngrdcol, dtype=float),
            l_implemented=True,
            l_ascending_grid=True,
            grid_type=2,
            deltaz=np.full(ngrdcol, float(np.median(np.diff(inputs["zt"]))), dtype=float),
            zm_init=np.asarray((inputs["zm"][0],), dtype=float),
            zm_top=np.asarray((inputs["zm"][-1],), dtype=float),
            momentum_heights=np.asfortranarray(inputs["zm"][None, :]),
            thermodynamic_heights=np.asfortranarray(inputs["zt"][None, :]),
            err_info=err_info,
        )
        params = api.init_clubb_params(ngrdcol, iunit=10, filename="")
        _nu, lmin, _mixture_limit = api.calc_derrived_params(
            gr=grid,
            ngrdcol=ngrdcol,
            grid_type=2,
            deltaz=np.full(ngrdcol, float(np.median(np.diff(inputs["zt"]))), dtype=float),
            clubb_params=params,
            nu_vert_res_dep=None,
            l_prescribed_avg_deltaz=bool(flags.l_prescribed_avg_deltaz),
        )
        pdf_params = api.init_pdf_params(nzt, ngrdcol)
        thv_ds = np.asfortranarray(inputs["thlm"][None, :])
        rt_path_2d = np.asfortranarray(rt_path[None, :])
        thl_path_2d = np.asfortranarray(thl_path[None, :])
        exner = np.asfortranarray(_exner(inputs["pressure"])[None, :])
        thvm = api.calculate_thvm(
            nzt,
            ngrdcol,
            thv_ds,
            np.asfortranarray(inputs["rtm"][None, :]),
            np.asfortranarray(inputs["rcm"][None, :]),
            exner,
            thv_ds,
        )
        err_info_out, lscale, lscale_up, lscale_down = api.calc_lscale_directly(
            ngrdcol=ngrdcol,
            nzm=nzm,
            nzt=nzt,
            gr=grid,
            l_implemented=True,
            p_in_pa=np.asfortranarray(inputs["pressure"][None, :]),
            exner=exner,
            rtm=rt_path_2d,
            thlm=thl_path_2d,
            thvm=np.asfortranarray(thvm),
            newmu=np.asarray((_parameter_value(api, params, "mu"),), dtype=float),
            rtp2_zt=np.asfortranarray(inputs["rtp2"][None, :]),
            thlp2_zt=np.asfortranarray(inputs["thlp2"][None, :]),
            rtpthlp_zt=np.asfortranarray(inputs["rtpthlp"][None, :]),
            pdf_params=pdf_params,
            em=np.asfortranarray(inputs["em"][None, :]),
            thv_ds_zt=thv_ds,
            lscale_max=np.asarray((1.0e5,), dtype=float),
            lmin=float(lmin),
            clubb_params=np.asfortranarray(params),
            saturation_formula=int(flags.saturation_formula),
            l_lscale_plume_centered=False,
            err_info=err_info,
        )
        errors = np.asarray(api.get_err_code(ngrdcol), dtype=int)
        try:
            api.cleanup_grid(grid)
            api.cleanup_err_info(err_info_out)
        except (AttributeError, RuntimeError):
            pass
    if np.any(errors != 0):
        raise SamLscaleError(f"calc_Lscale_directly returned CLUBB error code(s) {errors.tolist()}.")
    values = (np.asarray(lscale[0], dtype=float), np.asarray(lscale_up[0], dtype=float), np.asarray(lscale_down[0], dtype=float))
    if not all(np.all(np.isfinite(value)) for value in values):
        raise SamLscaleError("calc_Lscale_directly returned a non-finite Lscale profile.")
    return SamLscaleProfile(
        source_file=profile.source_file,
        elapsed_seconds=profile.elapsed_seconds,
        elapsed_minutes=profile.elapsed_minutes,
        height_m=np.asarray(profile.height_m, dtype=float),
        lscale_m=values[0],
        lscale_up_m=values[1],
        lscale_down_m=values[2],
        tke_proxy_m2_s2=np.asarray(inputs["em"], dtype=float),
        lmin_m=float(lmin),
        newmu_m_inv=_parameter_value(api, params, "mu"),
        path_kind=path_kind,
        approximation=approximation,
    )


@lru_cache(maxsize=8)
def compute_sam_lscale_profile(run_dir: str, elapsed_seconds: int) -> SamLscaleProfile:
    """Compute the ordinary mean-state direct Lscale profile for raw SAM."""

    profile = load_profile_moments(str(Path(run_dir).expanduser()), int(elapsed_seconds))
    return _compute_native_lscale_profile(
        profile,
        rt_path=np.asarray(profile.mean[:, 1], dtype=float),
        thl_path=np.asarray(profile.mean[:, 2], dtype=float),
        path_kind="mean",
        approximation=(
            "Resolved SAM thermodynamic means/covariances; em = 1.5 w'^2 "
            "because this micro output has no horizontal velocity components."
        ),
    )


def compute_sam_plume_lscale_profile(
    run_dir: str,
    elapsed_seconds: int,
    *,
    plume_rt: np.ndarray,
    plume_thl: np.ndarray,
) -> SamLscaleProfile:
    """Re-run native direct Lscale along a diagnosed G1 thermodynamic path.

    This deliberately does not reimplement CLUBB's entrainment or CAPE
    integration in Python.  The supplied G1 center curve replaces only the
    r_t/theta_l path passed to ``calc_Lscale_directly``; environmental virtual
    temperature, pressure, TKE, and all other inputs remain the raw-SAM
    grid-mean profiles.  It is a laboratory coupling, not the disabled
    production ``l_Lscale_plume_centered`` branch.
    """

    profile = load_profile_moments(str(Path(run_dir).expanduser()), int(elapsed_seconds))
    return _compute_native_lscale_profile(
        profile,
        rt_path=plume_rt,
        thl_path=plume_thl,
        path_kind="g1",
        approximation=(
            "Native direct-Lscale entrainment/CAPE calculation along the current "
            "PDF-10 G1 r_t/theta_l center path; environmental theta_v remains the "
            "resolved SAM mean and em = 1.5 w'^2 because U/V are unavailable."
        ),
    )


__all__ = [
    "SamLscaleError",
    "SamLscaleProfile",
    "compute_sam_lscale_profile",
    "compute_sam_plume_lscale_profile",
]
