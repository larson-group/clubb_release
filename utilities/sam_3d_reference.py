"""Small raw-SAM reader shared by active Dash diagnostics.

Experimental PDF fitting has its own archived support code.  Active
diagnostics only need immutable horizontal planes and CLUBB's linearized cloud
transform, so those dependencies live here as a compact shared reader.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from pathlib import Path
import re
from typing import Any

import netCDF4
import numpy as np
from scipy.special import ndtr

from utilities.les_chi_moments import derive_chi_moments


SHARED_3D_SAM_ROOT = Path(
    "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/"
    "JULY_2017_3D_RECREATIONS"
)
DEFAULT_SAM_RUN = SHARED_3D_SAM_ROOT / "arm_3d"
DEFAULT_BOMEX_SAM_RUN = SHARED_3D_SAM_ROOT / "bomex_3d"
_STEP_PATTERN = re.compile(r"_(\d{10})_micro\.nc$")


@dataclass
class RunInventory:
    run_dir: str
    steps_seconds: list[int]
    heights_m: list[float]
    first_minute: float
    last_minute: float
    snapshot_count: int


@dataclass
class SnapshotData:
    run_dir: str
    source_file: str
    elapsed_seconds: int
    elapsed_minutes: float
    absolute_time_days: float
    height_m: float
    pressure_pa: float
    sample_count: int
    samples: np.ndarray
    chi_samples: np.ndarray
    rc_samples: np.ndarray
    mean: np.ndarray
    covariance: np.ndarray
    wp3: float
    rtp3: float
    thlp3: float
    targets: np.ndarray
    cloud_fraction: float
    rcm: float
    w_in_cloud: float
    rcm_in_cloud: float
    chi_mean: float
    chi_variance: float
    chi_coeff_rt: float
    chi_coeff_thl: float
    cloud_cell_count: int
    half_transport_cell_count: int
    half_transport_domain_fraction: float
    half_transport_cloud_fraction: float
    tail_placement: dict[str, Any] | None = None


@dataclass
class ProfileMoments:
    """Raw-SAM grid-mean profiles needed by column diagnostics.

    The immutable SAM micro files carry the resolved ``w``, ``r_t``,
    ``theta_l``, and cloud-water samples but not CLUBB's diagnosed fields.
    This compact profile bundle is intentionally limited to moments that can
    be calculated directly from those resolved samples. ``mixed_third`` is
    ordered as ``<w'^2 r_t'>``, ``<w' r_t'^2>``, and
    ``<w' r_t' theta_l'>`` for Notes-only conditional-structure experiments.
    """

    run_dir: str
    source_file: str
    elapsed_seconds: int
    elapsed_minutes: float
    height_m: np.ndarray
    pressure_pa: np.ndarray
    mean: np.ndarray
    covariance: np.ndarray
    third: np.ndarray
    mixed_third: np.ndarray
    mean_rc: np.ndarray


@lru_cache(maxsize=32)
def _snapshot_files(run_dir: Path) -> list[tuple[int, Path]]:
    """Return physical elapsed seconds and files using NetCDF time metadata."""
    output_dir = Path(run_dir).expanduser().resolve() / "OUT_3D"
    timestep_files: list[tuple[int, Path]] = []
    for path in output_dir.glob("*_micro.nc"):
        match = _STEP_PATTERN.search(path.name)
        if match:
            timestep_files.append((int(match.group(1)), path))
    timestep_files.sort()
    if len(timestep_files) < 2:
        return timestep_files

    first_step, first_path = timestep_files[0]
    last_step, last_path = timestep_files[-1]
    step_span = last_step - first_step
    if step_span <= 0:
        return timestep_files
    try:
        first_time = _snapshot_time_seconds(first_path)
        last_time = _snapshot_time_seconds(last_path)
        seconds_per_step = (last_time - first_time) / step_span
    except (KeyError, OSError, ValueError):
        return timestep_files
    if not math.isfinite(seconds_per_step) or seconds_per_step <= 0.0:
        return timestep_files
    return [
        (int(round(timestep * seconds_per_step)), path)
        for timestep, path in timestep_files
    ]


def _snapshot_time_seconds(path: Path) -> float:
    """Read one SAM snapshot time coordinate in physical seconds."""
    with netCDF4.Dataset(path) as dataset:
        time = dataset.variables["time"]
        value = float(np.asarray(time[:], dtype=float).reshape(-1)[0])
        units = str(getattr(time, "units", "") or "").strip().lower()
    unit = units.split()[0] if units else ""
    factors = {
        "day": 86400.0,
        "days": 86400.0,
        "hour": 3600.0,
        "hours": 3600.0,
        "minute": 60.0,
        "minutes": 60.0,
        "second": 1.0,
        "seconds": 1.0,
    }
    if unit not in factors:
        raise ValueError(f"Unsupported SAM snapshot time units: {units or 'missing'}")
    return value * factors[unit]


def inventory_run(run_dir: Path = DEFAULT_SAM_RUN) -> RunInventory:
    """Return available times and heights without opening every snapshot."""
    pairs = _snapshot_files(Path(run_dir))
    if not pairs:
        raise FileNotFoundError(
            f"No SAM *_micro.nc files found under {Path(run_dir) / 'OUT_3D'}"
        )
    with netCDF4.Dataset(pairs[0][1]) as dataset:
        heights = np.asarray(dataset.variables["z"][:], dtype=float)
    steps = [step for step, _path in pairs]
    return RunInventory(
        run_dir=str(Path(run_dir).expanduser().resolve()),
        steps_seconds=steps,
        heights_m=heights.tolist(),
        first_minute=steps[0] / 60.0,
        last_minute=steps[-1] / 60.0,
        snapshot_count=len(steps),
    )


def _nearest_snapshot(run_dir: Path, elapsed_seconds: int) -> tuple[int, Path]:
    pairs = _snapshot_files(run_dir)
    if not pairs:
        raise FileNotFoundError(
            f"No SAM 3-D snapshots found in {Path(run_dir) / 'OUT_3D'}"
        )
    return min(pairs, key=lambda pair: abs(pair[0] - int(elapsed_seconds)))


@lru_cache(maxsize=32)
def load_snapshot_pressure_profile(
    run_dir: str,
    elapsed_seconds: int,
) -> tuple[int, np.ndarray, np.ndarray]:
    """Return the nearest snapshot's heights and pressure profile in Pa."""

    step, path = _nearest_snapshot(Path(run_dir).expanduser(), elapsed_seconds)
    with netCDF4.Dataset(path) as dataset:
        return (
            step,
            np.asarray(dataset.variables["z"][:], dtype=float),
            np.asarray(dataset.variables["p"][:], dtype=float) * 100.0,
        )


@lru_cache(maxsize=8)
def load_profile_moments(run_dir: str, elapsed_seconds: int) -> ProfileMoments:
    """Calculate resolved raw-SAM moments for every height in one snapshot.

    This is deliberately cached by immutable source snapshot.  Interactive
    consumers can therefore move a height slider without reopening every
    horizontal plane or redoing the full-column reductions.
    """

    step, path = _nearest_snapshot(Path(run_dir).expanduser(), elapsed_seconds)
    with netCDF4.Dataset(path) as dataset:
        heights = np.asarray(dataset.variables["z"][:], dtype=float)
        pressure_pa = np.asarray(dataset.variables["p"][:], dtype=float) * 100.0
        nz = int(heights.size)
        mean = np.empty((nz, 3), dtype=float)
        covariance = np.empty((nz, 3, 3), dtype=float)
        third = np.empty((nz, 3), dtype=float)
        mixed_third = np.empty((nz, 3), dtype=float)
        mean_rc = np.empty(nz, dtype=float)
        for level in range(nz):
            samples = np.column_stack(
                tuple(
                    np.asarray(dataset.variables[name][0, level], dtype=float).reshape(-1)
                    for name in ("W", "RT", "THL")
                )
            )
            level_mean = np.mean(samples, axis=0)
            mean[level] = level_mean
            centered = samples - level_mean
            covariance[level] = _safe_covariance(samples)
            third[level] = np.mean(centered**3, axis=0)
            mixed_third[level] = np.asarray(
                (
                    np.mean(centered[:, 0] ** 2 * centered[:, 1]),
                    np.mean(centered[:, 0] * centered[:, 1] ** 2),
                    np.mean(centered[:, 0] * centered[:, 1] * centered[:, 2]),
                ),
                dtype=float,
            )
            mean_rc[level] = float(np.mean(dataset.variables["RC"][0, level]))

    return ProfileMoments(
        run_dir=str(Path(run_dir).expanduser().resolve()),
        source_file=str(path),
        elapsed_seconds=int(step),
        elapsed_minutes=float(step) / 60.0,
        height_m=heights,
        pressure_pa=pressure_pa,
        mean=mean,
        covariance=covariance,
        third=third,
        mixed_third=mixed_third,
        mean_rc=mean_rc,
    )


def _safe_covariance(samples: np.ndarray) -> np.ndarray:
    covariance = np.cov(np.asarray(samples, dtype=float), rowvar=False, bias=True)
    covariance = 0.5 * (covariance + covariance.T)
    values, vectors = np.linalg.eigh(covariance)
    floor = max(float(np.max(np.abs(values))), 1.0e-20) * 1.0e-11
    return vectors @ np.diag(np.maximum(values, floor)) @ vectors.T


def load_snapshot(
    run_dir: Path,
    elapsed_seconds: int,
    height_m: float,
) -> SnapshotData:
    """Load one raw SAM plane and calculate CLUBB-style central moments."""
    step, path = _nearest_snapshot(Path(run_dir).expanduser(), elapsed_seconds)
    with netCDF4.Dataset(path) as dataset:
        z = np.asarray(dataset.variables["z"][:], dtype=float)
        level = int(np.argmin(np.abs(z - float(height_m))))
        arrays = {
            name: np.asarray(dataset.variables[name][0, level], dtype=float).reshape(-1)
            for name in ("W", "RT", "THL", "CHI", "RC")
        }
        pressure_pa = float(dataset.variables["p"][level]) * 100.0
        absolute_time_days = float(dataset.variables["time"][0])

    samples = np.column_stack((arrays["W"], arrays["RT"], arrays["THL"]))
    mean = np.mean(samples, axis=0)
    centered = samples - mean
    covariance = _safe_covariance(samples)
    rc_centered = arrays["RC"] - float(np.mean(arrays["RC"]))
    targets = np.array(
        (
            np.mean(centered[:, 0] * rc_centered),
            np.mean(centered[:, 0] ** 2 * rc_centered),
            np.mean(centered[:, 1] * rc_centered),
            np.mean(centered[:, 2] * rc_centered),
        ),
        dtype=float,
    )
    cloudy = arrays["RC"] > 0.0
    cloud_count = int(np.count_nonzero(cloudy))
    cloud_fraction = float(np.mean(cloudy))
    rcm = float(np.mean(arrays["RC"]))
    w_in_cloud = float(np.mean(arrays["W"][cloudy])) if cloud_count else math.nan
    rcm_in_cloud = float(np.mean(arrays["RC"][cloudy])) if cloud_count else math.nan
    chi = derive_chi_moments(
        mean[1],
        mean[2],
        pressure_pa,
        covariance[1, 1],
        covariance[2, 2],
        covariance[1, 2],
        covariance[0, 1],
        covariance[0, 2],
        covariance[0, 0],
    )
    local_transport = centered[:, 0] * arrays["RC"]
    order = np.argsort(local_transport)[::-1]
    net = float(np.sum(local_transport))
    half_count = (
        int(np.searchsorted(np.cumsum(local_transport[order]), 0.5 * net) + 1)
        if net > 0.0
        else 0
    )
    return SnapshotData(
        run_dir=str(Path(run_dir).expanduser().resolve()),
        source_file=str(path),
        elapsed_seconds=step,
        elapsed_minutes=step / 60.0,
        absolute_time_days=absolute_time_days,
        height_m=float(z[level]),
        pressure_pa=pressure_pa,
        sample_count=int(samples.shape[0]),
        samples=samples,
        chi_samples=arrays["CHI"],
        rc_samples=arrays["RC"],
        mean=mean,
        covariance=covariance,
        wp3=float(np.mean(centered[:, 0] ** 3)),
        rtp3=float(np.mean(centered[:, 1] ** 3)),
        thlp3=float(np.mean(centered[:, 2] ** 3)),
        targets=targets,
        cloud_fraction=cloud_fraction,
        rcm=rcm,
        w_in_cloud=w_in_cloud,
        rcm_in_cloud=rcm_in_cloud,
        chi_mean=float(chi["mean_chi"]),
        chi_variance=float(chi["var_chi"]),
        chi_coeff_rt=float(chi["coef_rt"]),
        chi_coeff_thl=float(chi["coef_thl"]),
        cloud_cell_count=cloud_count,
        half_transport_cell_count=half_count,
        half_transport_domain_fraction=half_count / samples.shape[0],
        half_transport_cloud_fraction=(
            half_count / cloud_count if cloud_count else math.nan
        ),
    )


def diagnose_cloud(
    weights: np.ndarray,
    means: np.ndarray,
    covariances: np.ndarray,
    snapshot: SnapshotData,
) -> dict[str, Any]:
    """Evaluate CLUBB's linearized, truncated-normal cloud integrals."""
    weights = np.asarray(weights, dtype=float)
    means = np.asarray(means, dtype=float)
    covariances = np.asarray(covariances, dtype=float)
    displacement = means - snapshot.mean[None, :]
    crt = snapshot.chi_coeff_rt
    cthl = snapshot.chi_coeff_thl
    chi_mean = snapshot.chi_mean + crt * displacement[:, 1] - cthl * displacement[:, 2]
    chi_variance = np.maximum(
        crt**2 * covariances[:, 1, 1]
        + cthl**2 * covariances[:, 2, 2]
        - 2.0 * crt * cthl * covariances[:, 1, 2],
        1.0e-30,
    )
    sigma_chi = np.sqrt(chi_variance)
    cov_w_chi = crt * covariances[:, 0, 1] - cthl * covariances[:, 0, 2]
    cov_rt_chi = crt * covariances[:, 1, 1] - cthl * covariances[:, 1, 2]
    cov_thl_chi = crt * covariances[:, 1, 2] - cthl * covariances[:, 2, 2]
    normalized = chi_mean / sigma_chi
    cloud_component = ndtr(normalized)
    normal_pdf = np.exp(-0.5 * np.clip(normalized, -40.0, 40.0) ** 2) / math.sqrt(
        2.0 * math.pi
    )
    rc_component = chi_mean * cloud_component + sigma_chi * normal_pdf
    cloud_fraction = float(np.sum(weights * cloud_component))
    rcm = float(np.sum(weights * rc_component))
    centered_rc = rc_component - rcm
    wprcp = np.sum(
        weights * (displacement[:, 0] * centered_rc + cov_w_chi * cloud_component)
    )
    wp2rcp = np.sum(
        weights
        * (
            (displacement[:, 0] ** 2 + covariances[:, 0, 0]) * centered_rc
            + 2.0 * displacement[:, 0] * cov_w_chi * cloud_component
            + cov_w_chi**2 / sigma_chi * normal_pdf
        )
    )
    rtprcp = np.sum(
        weights * (displacement[:, 1] * centered_rc + cov_rt_chi * cloud_component)
    )
    thlprcp = np.sum(
        weights * (displacement[:, 2] * centered_rc + cov_thl_chi * cloud_component)
    )
    w_cloud_numerator = np.sum(
        weights * (means[:, 0] * cloud_component + cov_w_chi / sigma_chi * normal_pdf)
    )
    return {
        "diagnostics": np.array((wprcp, wp2rcp, rtprcp, thlprcp), dtype=float),
        "cloud_fraction": cloud_fraction,
        "rcm": rcm,
        "w_in_cloud": w_cloud_numerator / max(cloud_fraction, 1.0e-12),
        "rcm_in_cloud": rcm / max(cloud_fraction, 1.0e-12),
        "chi_mean_component": chi_mean,
        "chi_variance_component": chi_variance,
        "cov_w_chi_component": cov_w_chi,
        "cloud_fraction_component": cloud_component,
    }
