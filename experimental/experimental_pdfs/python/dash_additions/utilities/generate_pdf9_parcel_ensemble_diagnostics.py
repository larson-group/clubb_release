#!/usr/bin/env python3
"""Build the compact cache used by the PDF-9 parcel-ensemble Misc lab.

The Fortran PDF-9 experiment writes a launch-by-destination ledger.  This
utility copies that immutable ledger into a compressed, browser-friendly cube
and reduces the matching raw SAM planes to evaluation-only centers, shapes,
and cloud diagnostics.  No SAM conditional statistic is supplied to CLUBB or
used to evolve a parcel.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Iterable

import netCDF4
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "output_pdf_bakeoff/pdf9_parcel_ensemble_lab"
SAM_ROOT = Path(
    "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/"
    "JULY_2017_3D_RECREATIONS"
)
SCHEMA_VERSION = 16
STATE_NAMES = ("w", "rt", "thl")
REFERENCE_NAMES = (
    "exact",
    "all_cloud",
    "cloud_up",
    "cloud_down",
    "clear_moist_up",
)


def _first_existing(paths: Iterable[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    return next(iter(paths))


DEFAULT_CASES = {
    "arm": (
        _first_existing(
            (
                REPO_ROOT
                / "output_pdf_bakeoff/pdf9_parcel_runs/weighted/arm/arm_stats.nc",
                REPO_ROOT
                / "output_pdf_bakeoff/pdf9_parcel_runs/recursive/arm/arm_stats.nc",
                REPO_ROOT
                / "output_pdf_bakeoff/pdf9_parcel_runs/arm/arm_stats.nc",
                Path("/tmp/pdf9_stepbatch/arm_stats.nc"),
                REPO_ROOT / "output/arm_stats.nc",
            )
        ),
        SAM_ROOT / "arm_3d",
    ),
    "bomex": (
        _first_existing(
            (
                REPO_ROOT
                / "output_pdf_bakeoff/pdf9_parcel_runs/weighted/bomex/bomex_stats.nc",
                REPO_ROOT
                / "output_pdf_bakeoff/pdf9_parcel_runs/recursive/bomex/bomex_stats.nc",
                Path("/tmp/pdf9_stepbatch_bomex/bomex_stats.nc"),
                REPO_ROOT / "output/bomex_stats.nc",
            )
        ),
        SAM_ROOT / "bomex_3d",
    ),
}


def _interpolate_time_height(
    values: np.ndarray,
    source_height: np.ndarray,
    destination_height: np.ndarray,
) -> np.ndarray:
    """Interpolate a time-height field without smearing across time."""
    result = np.empty((values.shape[0], len(destination_height)), np.float32)
    for time_index, profile in enumerate(values):
        result[time_index] = np.interp(
            destination_height, source_height, profile
        ).astype(np.float32)
    return result


def _weighted_moments(samples: np.ndarray, weights: np.ndarray):
    weights = np.asarray(weights, dtype=float).reshape(-1)
    samples = np.asarray(samples, dtype=float).reshape(-1, 3)
    finite = np.isfinite(weights) & (weights > 0.0) & np.all(
        np.isfinite(samples), axis=1
    )
    total = float(np.sum(weights[finite]))
    if total <= 1.0e-24:
        return False, np.full(3, np.nan), np.full((3, 3), np.nan), 0.0
    normalized = weights[finite] / total
    selected = samples[finite]
    mean = np.sum(normalized[:, None] * selected, axis=0)
    anomaly = selected - mean
    covariance = np.einsum("n,ni,nj->ij", normalized, anomaly, anomaly)
    return True, mean, covariance, total / float(len(weights))


def _snapshot_inventory(run_dir: Path) -> list[tuple[int, Path]]:
    pattern = re.compile(r"_(\d+)_micro\.nc$")
    result = []
    for path in (run_dir / "OUT_3D").glob("*_micro.nc"):
        match = pattern.search(path.name)
        if match:
            result.append((int(match.group(1)), path))
    if not result:
        raise FileNotFoundError(f"No raw SAM snapshots under {run_dir / 'OUT_3D'}")
    return sorted(result)


def _nearest_file(inventory: list[tuple[int, Path]], seconds: float):
    return min(inventory, key=lambda item: abs(item[0] - seconds))


def _read_stats(path: Path, raw_max_height: float) -> dict[str, np.ndarray]:
    required = (
        "pdf9_up_parcel_tke",
        "pdf9_up_parcel_w",
        "pdf9_up_parcel_rt",
        "pdf9_up_parcel_thl",
        "pdf9_up_parcel_buoyancy",
        "pdf9_up_parcel_status",
        "pdf9_down_parcel_tke",
        "pdf9_down_parcel_w",
        "pdf9_down_parcel_rt",
        "pdf9_down_parcel_thl",
        "pdf9_down_parcel_buoyancy",
        "pdf9_down_parcel_status",
        "pdf9_lscale_up",
        "pdf9_lscale_down",
        "pdf9_candidate_valid_up",
        "pdf9_candidate_launch_from_g3_up",
        "pdf9_candidate_launch_from_g3_down",
        "pdf9_candidate_crossing_count_up",
        "pdf9_candidate_weighted_support_up",
        "pdf9_candidate_donor_distance_up",
        "pdf9_candidate_w_up",
        "pdf9_candidate_rt_up",
        "pdf9_candidate_thl_up",
        "pdf9_candidate_var_w_up",
        "pdf9_candidate_var_rt_up",
        "pdf9_candidate_var_thl_up",
        "pdf9_candidate_covar_w_rt_up",
        "pdf9_candidate_covar_w_thl_up",
        "pdf9_candidate_covar_rt_thl_up",
        "pdf9_candidate_corr_w_rt_up",
        "pdf9_candidate_corr_w_thl_up",
        "pdf9_candidate_corr_rt_thl_up",
        "pdf9_candidate_valid_down",
        "pdf9_candidate_crossing_count_down",
        "pdf9_candidate_weighted_support_down",
        "pdf9_candidate_donor_distance_down",
        "pdf9_candidate_w_down",
        "pdf9_candidate_w_down_uncapped",
        "pdf9_candidate_down_cap_fraction",
        "pdf9_candidate_destination_sigma_w",
        "pdf9_candidate_rt_down",
        "pdf9_candidate_thl_down",
        "pdf9_candidate_var_w_down",
        "pdf9_candidate_var_rt_down",
        "pdf9_candidate_var_thl_down",
        "pdf9_candidate_covar_w_rt_down",
        "pdf9_candidate_covar_w_thl_down",
        "pdf9_candidate_covar_rt_thl_down",
        "pdf9_candidate_corr_w_rt_down",
        "pdf9_candidate_corr_w_thl_down",
        "pdf9_candidate_corr_rt_thl_down",
        "pdf9_candidate_branch_prob_up",
        "pdf9_candidate_branch_prob_down",
        "pdf9_candidate_valid_combined",
        "pdf9_candidate_crossing_count_combined",
        "pdf9_candidate_weighted_support_combined",
        "pdf9_candidate_donor_distance_combined",
        "pdf9_candidate_w_combined",
        "pdf9_candidate_rt_combined",
        "pdf9_candidate_thl_combined",
        "pdf9_candidate_var_w_combined",
        "pdf9_candidate_var_rt_combined",
        "pdf9_candidate_var_thl_combined",
        "pdf9_candidate_covar_w_rt_combined",
        "pdf9_candidate_covar_w_thl_combined",
        "pdf9_candidate_covar_rt_thl_combined",
        "pdf9_candidate_corr_w_rt_combined",
        "pdf9_candidate_corr_w_thl_combined",
        "pdf9_candidate_corr_rt_thl_combined",
        "mixt_frac",
        "w_1",
        "w_2",
        "rt_1",
        "rt_2",
        "varnce_w_1",
        "varnce_w_2",
        "varnce_rt_1",
        "varnce_rt_2",
        "corr_w_rt_1",
        "corr_w_rt_2",
        "chi_1",
        "chi_2",
        "stdev_chi_1",
        "stdev_chi_2",
        "corr_w_chi_1",
        "corr_w_chi_2",
        "rc_1",
        "rc_2",
        "mixt_frac_3",
        "w_3",
        "rt_3",
        "thl_3",
        "varnce_w_3",
        "varnce_rt_3",
        "varnce_thl_3",
        "corr_w_rt_3",
        "corr_w_thl_3",
        "corr_rt_thl_3",
        "chi_3",
        "stdev_chi_3",
        "corr_w_chi_3",
        "rc_3",
        "rtm",
        "thlm",
        "pdf9_entrain_g3_weight",
        "pdf9_entrain_rt_up",
        "pdf9_entrain_rt_down",
        "pdf9_entrain_thl_up",
        "pdf9_entrain_thl_down",
    )
    with netCDF4.Dataset(path) as dataset:
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            raise ValueError(
                f"{path} does not contain the PDF-9 parcel ledger: {', '.join(missing)}"
            )
        destination = np.asarray(dataset.variables["pdf9_destination_zt"][:], float)
        launch = np.asarray(dataset.variables["pdf9_launch_zt"][:], float)
        keep = np.where(destination <= raw_max_height + 1.0e-6)[0]
        if not keep.size:
            raise ValueError("No PDF-9 destination levels overlap raw SAM")
        data = {
            "time_seconds": np.asarray(dataset.variables["time"][:], float),
            "height_m": destination[keep],
            "launch_height_m": launch,
        }
        for short, source in (
            ("parcel_tke", "pdf9_up_parcel_tke"),
            ("parcel_w", "pdf9_up_parcel_w"),
            ("parcel_rt", "pdf9_up_parcel_rt"),
            ("parcel_thl", "pdf9_up_parcel_thl"),
            ("parcel_buoyancy", "pdf9_up_parcel_buoyancy"),
        ):
            data[short] = np.asarray(dataset.variables[source][:, keep, :, 0], np.float32)
        data["parcel_status"] = np.asarray(
            dataset.variables["pdf9_up_parcel_status"][:, keep, :, 0], np.int8
        )
        for short, source in (
            ("down_parcel_tke", "pdf9_down_parcel_tke"),
            ("down_parcel_w", "pdf9_down_parcel_w"),
            ("down_parcel_rt", "pdf9_down_parcel_rt"),
            ("down_parcel_thl", "pdf9_down_parcel_thl"),
            ("down_parcel_buoyancy", "pdf9_down_parcel_buoyancy"),
        ):
            data[short] = np.asarray(
                dataset.variables[source][:, keep, :, 0], np.float32
            )
        data["down_parcel_status"] = np.asarray(
            dataset.variables["pdf9_down_parcel_status"][:, keep, :, 0], np.int8
        )
        for short, source in (
            ("lscale_up", "pdf9_lscale_up"),
            ("lscale_down", "pdf9_lscale_down"),
            ("candidate_valid_up", "pdf9_candidate_valid_up"),
            (
                "candidate_launch_from_g3_up",
                "pdf9_candidate_launch_from_g3_up",
            ),
            (
                "candidate_launch_from_g3_down",
                "pdf9_candidate_launch_from_g3_down",
            ),
            (
                "candidate_crossing_count_up",
                "pdf9_candidate_crossing_count_up",
            ),
            (
                "candidate_weighted_support_up",
                "pdf9_candidate_weighted_support_up",
            ),
            (
                "candidate_donor_distance_up",
                "pdf9_candidate_donor_distance_up",
            ),
            ("candidate_w_up", "pdf9_candidate_w_up"),
            ("candidate_rt_up", "pdf9_candidate_rt_up"),
            ("candidate_thl_up", "pdf9_candidate_thl_up"),
            ("candidate_var_w_up", "pdf9_candidate_var_w_up"),
            ("candidate_var_rt_up", "pdf9_candidate_var_rt_up"),
            ("candidate_var_thl_up", "pdf9_candidate_var_thl_up"),
            ("candidate_covar_w_rt_up", "pdf9_candidate_covar_w_rt_up"),
            ("candidate_covar_w_thl_up", "pdf9_candidate_covar_w_thl_up"),
            ("candidate_covar_rt_thl_up", "pdf9_candidate_covar_rt_thl_up"),
            ("candidate_corr_w_rt_up", "pdf9_candidate_corr_w_rt_up"),
            ("candidate_corr_w_thl_up", "pdf9_candidate_corr_w_thl_up"),
            ("candidate_corr_rt_thl_up", "pdf9_candidate_corr_rt_thl_up"),
            ("candidate_valid_down", "pdf9_candidate_valid_down"),
            (
                "candidate_crossing_count_down",
                "pdf9_candidate_crossing_count_down",
            ),
            (
                "candidate_weighted_support_down",
                "pdf9_candidate_weighted_support_down",
            ),
            (
                "candidate_donor_distance_down",
                "pdf9_candidate_donor_distance_down",
            ),
            ("candidate_w_down", "pdf9_candidate_w_down"),
            (
                "candidate_w_down_uncapped",
                "pdf9_candidate_w_down_uncapped",
            ),
            (
                "candidate_down_cap_fraction",
                "pdf9_candidate_down_cap_fraction",
            ),
            (
                "candidate_destination_sigma_w",
                "pdf9_candidate_destination_sigma_w",
            ),
            ("candidate_rt_down", "pdf9_candidate_rt_down"),
            ("candidate_thl_down", "pdf9_candidate_thl_down"),
            ("candidate_var_w_down", "pdf9_candidate_var_w_down"),
            ("candidate_var_rt_down", "pdf9_candidate_var_rt_down"),
            ("candidate_var_thl_down", "pdf9_candidate_var_thl_down"),
            ("candidate_covar_w_rt_down", "pdf9_candidate_covar_w_rt_down"),
            ("candidate_covar_w_thl_down", "pdf9_candidate_covar_w_thl_down"),
            (
                "candidate_covar_rt_thl_down",
                "pdf9_candidate_covar_rt_thl_down",
            ),
            ("candidate_corr_w_rt_down", "pdf9_candidate_corr_w_rt_down"),
            ("candidate_corr_w_thl_down", "pdf9_candidate_corr_w_thl_down"),
            ("candidate_corr_rt_thl_down", "pdf9_candidate_corr_rt_thl_down"),
            ("candidate_branch_prob_up", "pdf9_candidate_branch_prob_up"),
            ("candidate_branch_prob_down", "pdf9_candidate_branch_prob_down"),
            ("candidate_valid_combined", "pdf9_candidate_valid_combined"),
            (
                "candidate_crossing_count_combined",
                "pdf9_candidate_crossing_count_combined",
            ),
            (
                "candidate_weighted_support_combined",
                "pdf9_candidate_weighted_support_combined",
            ),
            (
                "candidate_donor_distance_combined",
                "pdf9_candidate_donor_distance_combined",
            ),
            ("candidate_w_combined", "pdf9_candidate_w_combined"),
            ("candidate_rt_combined", "pdf9_candidate_rt_combined"),
            ("candidate_thl_combined", "pdf9_candidate_thl_combined"),
            ("candidate_var_w_combined", "pdf9_candidate_var_w_combined"),
            ("candidate_var_rt_combined", "pdf9_candidate_var_rt_combined"),
            ("candidate_var_thl_combined", "pdf9_candidate_var_thl_combined"),
            (
                "candidate_covar_w_rt_combined",
                "pdf9_candidate_covar_w_rt_combined",
            ),
            (
                "candidate_covar_w_thl_combined",
                "pdf9_candidate_covar_w_thl_combined",
            ),
            (
                "candidate_covar_rt_thl_combined",
                "pdf9_candidate_covar_rt_thl_combined",
            ),
            ("candidate_corr_w_rt_combined", "pdf9_candidate_corr_w_rt_combined"),
            (
                "candidate_corr_w_thl_combined",
                "pdf9_candidate_corr_w_thl_combined",
            ),
            (
                "candidate_corr_rt_thl_combined",
                "pdf9_candidate_corr_rt_thl_combined",
            ),
            ("g12_mixt_frac", "mixt_frac"),
            ("g1_w", "w_1"),
            ("g2_w", "w_2"),
            ("g1_rt", "rt_1"),
            ("g2_rt", "rt_2"),
            ("g1_var_w", "varnce_w_1"),
            ("g2_var_w", "varnce_w_2"),
            ("g1_var_rt", "varnce_rt_1"),
            ("g2_var_rt", "varnce_rt_2"),
            ("g1_corr_w_rt", "corr_w_rt_1"),
            ("g2_corr_w_rt", "corr_w_rt_2"),
            ("g1_chi", "chi_1"),
            ("g2_chi", "chi_2"),
            ("g1_stdev_chi", "stdev_chi_1"),
            ("g2_stdev_chi", "stdev_chi_2"),
            ("g1_corr_w_chi", "corr_w_chi_1"),
            ("g2_corr_w_chi", "corr_w_chi_2"),
            ("g1_rc", "rc_1"),
            ("g2_rc", "rc_2"),
            ("g3_weight", "mixt_frac_3"),
            ("g3_w", "w_3"),
            ("g3_rt", "rt_3"),
            ("g3_thl", "thl_3"),
            ("g3_var_w", "varnce_w_3"),
            ("g3_var_rt", "varnce_rt_3"),
            ("g3_var_thl", "varnce_thl_3"),
            ("g3_corr_w_rt", "corr_w_rt_3"),
            ("g3_corr_w_thl", "corr_w_thl_3"),
            ("g3_corr_rt_thl", "corr_rt_thl_3"),
            ("g3_chi", "chi_3"),
            ("g3_stdev_chi", "stdev_chi_3"),
            ("g3_corr_w_chi", "corr_w_chi_3"),
            ("g3_rc", "rc_3"),
            ("mean_rt", "rtm"),
            ("mean_thl", "thlm"),
            ("entrain_g3_weight", "pdf9_entrain_g3_weight"),
            ("entrain_rt_up", "pdf9_entrain_rt_up"),
            ("entrain_rt_down", "pdf9_entrain_rt_down"),
            ("entrain_thl_up", "pdf9_entrain_thl_up"),
            ("entrain_thl_down", "pdf9_entrain_thl_down"),
        ):
            values = dataset.variables[source][:, keep, 0]
            fill = (
                0.0
                if short in ("g3_weight", "entrain_g3_weight")
                or short.startswith("g3_corr")
                or short.startswith("g1_corr")
                or short.startswith("g2_corr")
                else np.nan
            )
            data[short] = np.asarray(np.ma.filled(values, fill), np.float32)
        data["mean_w"] = np.zeros_like(data["mean_rt"], dtype=np.float32)

        # The parcel ledger stores deterministic crossing centers, not an
        # internal covariance for each parcel.  Preserve the local CLUBB
        # covariance at every destination so the report can draw a deliberately
        # small, clearly labelled footprint around each center and can show the
        # actual grid-mean r_t +/- 2 sigma envelope.
        zm = np.asarray(dataset.variables["zm"][:], float)
        destination_kept = destination[keep]
        for short, source in (
            ("local_var_w", "wp2"),
            ("local_var_rt", "rtp2"),
            ("local_cov_w_rt", "wprtp"),
        ):
            values = np.asarray(dataset.variables[source][:, :, 0], float)
            data[short] = _interpolate_time_height(
                values, zm, destination_kept
            )
    return data


def _reduce_sam(
    run_dir: Path,
    requested_time: np.ndarray,
    requested_height: np.ndarray,
) -> dict[str, np.ndarray]:
    inventory = _snapshot_inventory(run_dir)
    with netCDF4.Dataset(inventory[0][1]) as source:
        raw_heights = np.asarray(source.variables["z"][:], float)
    level_indices = np.asarray(
        [int(np.argmin(np.abs(raw_heights - height))) for height in requested_height]
    )
    shape = (len(requested_time), len(requested_height))
    means = np.full(shape + (3,), np.nan, np.float32)
    covariances = np.full(shape + (3, 3), np.nan, np.float32)
    reference_mean = np.full(
        shape + (len(REFERENCE_NAMES), 3), np.nan, np.float32
    )
    reference_covariance = np.full(
        shape + (len(REFERENCE_NAMES), 3, 3), np.nan, np.float32
    )
    reference_strength = np.zeros(shape + (len(REFERENCE_NAMES),), np.float32)
    reference_available = np.zeros(shape + (len(REFERENCE_NAMES),), np.int8)
    targets = np.full(shape + (4,), np.nan, np.float32)
    cloud_fraction = np.zeros(shape, np.float32)
    rcm = np.zeros(shape, np.float32)
    matched_seconds = np.zeros(len(requested_time), np.float64)

    for time_index, seconds in enumerate(requested_time):
        if time_index % 100 == 0:
            print(
                f"  raw SAM plane {time_index + 1}/{len(requested_time)}",
                flush=True,
            )
        raw_seconds, path = _nearest_file(inventory, float(seconds))
        matched_seconds[time_index] = raw_seconds
        with netCDF4.Dataset(path) as source:
            # One slab read per variable is much faster than thousands of
            # height-by-height netCDF reads during a full ARM reduction.
            slabs = {
                name: np.asarray(source.variables[name][0, level_indices], float)
                .reshape(len(level_indices), -1)
                for name in ("W", "RT", "THL", "RC")
            }
            for height_index, _level in enumerate(level_indices):
                arrays = {name: values[height_index] for name, values in slabs.items()}
                samples = np.column_stack((arrays["W"], arrays["RT"], arrays["THL"]))
                mean = np.mean(samples, axis=0)
                anomaly = samples - mean
                rc = arrays["RC"]
                w_prime = anomaly[:, 0]
                rt_prime = anomaly[:, 1]
                weights = (
                    np.ones(len(samples)),
                    rc,
                    rc * np.maximum(w_prime, 0.0),
                    rc * np.maximum(-w_prime, 0.0),
                    (rc <= 0.0)
                    * np.maximum(w_prime, 0.0)
                    * np.maximum(rt_prime, 0.0),
                )
                means[time_index, height_index] = mean
                covariances[time_index, height_index] = np.cov(
                    samples, rowvar=False, bias=True
                )
                for reference_index, weight in enumerate(weights):
                    available, center, covariance, strength = _weighted_moments(
                        samples, weight
                    )
                    reference_available[
                        time_index, height_index, reference_index
                    ] = available
                    reference_mean[
                        time_index, height_index, reference_index
                    ] = center
                    reference_covariance[
                        time_index, height_index, reference_index
                    ] = covariance
                    reference_strength[
                        time_index, height_index, reference_index
                    ] = strength
                rc_prime = rc - np.mean(rc)
                targets[time_index, height_index] = (
                    np.mean(anomaly[:, 0] * rc_prime),
                    np.mean(anomaly[:, 0] ** 2 * rc_prime),
                    np.mean(anomaly[:, 1] * rc_prime),
                    np.mean(anomaly[:, 2] * rc_prime),
                )
                cloud_fraction[time_index, height_index] = np.mean(rc > 0.0)
                rcm[time_index, height_index] = np.mean(rc)
    return {
        "sam_time_seconds": matched_seconds,
        "sam_height_m": raw_heights[level_indices].astype(np.float32),
        "sam_mean": means,
        "sam_covariance": covariances,
        "sam_reference_mean": reference_mean,
        "sam_reference_covariance": reference_covariance,
        "sam_reference_strength": reference_strength,
        "sam_reference_available": reference_available,
        "sam_targets": targets,
        "sam_cloud_fraction": cloud_fraction,
        "sam_rcm": rcm,
    }


def _cached_sam_reduction(
    path: Path,
    requested_time: np.ndarray,
    requested_height: np.ndarray,
) -> dict[str, np.ndarray] | None:
    """Reuse the immutable raw-SAM reduction when its coordinates match."""
    required = (
        "sam_time_seconds",
        "sam_height_m",
        "sam_mean",
        "sam_covariance",
        "sam_reference_mean",
        "sam_reference_covariance",
        "sam_reference_strength",
        "sam_reference_available",
        "sam_targets",
        "sam_cloud_fraction",
        "sam_rcm",
    )
    if not path.exists():
        return None
    try:
        with np.load(path, allow_pickle=False) as source:
            if any(name not in source for name in required):
                return None
            if not np.allclose(
                source["time_seconds"], requested_time, rtol=0.0, atol=1.0e-9
            ):
                return None
            if not np.allclose(
                source["height_m"], requested_height, rtol=0.0, atol=1.0e-6
            ):
                return None
            return {name: np.asarray(source[name]).copy() for name in required}
    except (KeyError, OSError, ValueError):
        return None


def build_case(name: str, stats_path: Path, sam_run: Path, output_dir: Path) -> dict:
    inventory = _snapshot_inventory(sam_run)
    with netCDF4.Dataset(inventory[0][1]) as source:
        raw_max_height = float(np.max(source.variables["z"][:]))
    cube = _read_stats(stats_path, raw_max_height)
    # CLUBB case time may include a case-specific start offset (ARM starts at
    # 41,400 s), whereas SAM snapshot filenames are elapsed LES seconds.  Both
    # streams are output after the first interval, so align their first saved
    # records and retain the absolute CLUBB clock separately for auditing.
    model_time = np.asarray(cube["time_seconds"], float).copy()
    sam_time = model_time - model_time[0] + float(inventory[0][0])
    cube["model_time_seconds"] = model_time
    cube["time_seconds"] = sam_time
    data_file = f"{name}.npz"
    data_path = output_dir / data_file
    sam_reduction = _cached_sam_reduction(
        data_path, sam_time, cube["height_m"]
    )
    if sam_reduction is None:
        sam_reduction = _reduce_sam(sam_run, sam_time, cube["height_m"])
    else:
            print(f"Reusing immutable SAM reduction from {data_path}", flush=True)
    cube.update(sam_reduction)
    temporary = data_path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **cube)
    temporary.replace(data_path)
    return {
        "name": name,
        "data_file": data_file,
        "stats_path": str(stats_path.resolve()),
        "sam_run_dir": str(sam_run.resolve()),
        "time_count": int(len(cube["time_seconds"])),
        "height_count": int(len(cube["height_m"])),
        "launch_count": int(len(cube["launch_height_m"])),
    }


def _parse_case(value: str):
    parts = value.split("::")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("Use NAME::STATS_PATH::SAM_RUN_DIR")
    return parts[0], Path(parts[1]).expanduser(), Path(parts[2]).expanduser()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        type=_parse_case,
        help="Repeatable NAME::STATS_PATH::SAM_RUN_DIR specification",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    cases = args.case or [
        (name, paths[0], paths[1]) for name, paths in DEFAULT_CASES.items()
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state_names": list(STATE_NAMES),
        "reference_names": list(REFERENCE_NAMES),
        "cases": [],
        "notes": {
            "sam_role": "evaluation only",
            "crossing_weight": "not a diagnosed PDF mixture fraction",
            "recursive_candidate": (
                "component-3 geometry is persisted for the next upward and "
                "downward parcel launches; the current candidate receives an "
                "analytic w-prime-r-t-prime transport weight"
            ),
            "downward_state": (
                "full signed-w launch-by-destination ledger; speed magnitude "
                "is available for diagnostic aggregation"
            ),
            "directional_candidates": (
                "upward and downward crossings are weighted by "
                "branch probability times exp(-mu times donor distance), "
                "stored separately, and pooled "
                "exactly into the next G3 candidate; directional launches use the "
                "positive- and negative-w conditional halves, capped at one "
                "local grid-mean vertical-velocity RMS"
            ),
            "downward_arrival_cap": (
                "downward trajectories and reach remain raw, but their velocity "
                "contribution to the next G3 mean and covariance is capped at "
                "minus the destination-local sqrt(wp2)"
            ),
            "entrainment_environment": (
                "parcel buoyancy entrainment blends the local grid mean "
                "directly with the directional G3 branch using the active "
                "G3 probability; launch and stopping rules are unchanged"
            ),
        },
    }
    for name, stats_path, sam_run in cases:
        if not stats_path.exists():
            print(f"Skipping {name}: no statistics file at {stats_path}")
            continue
        print(f"Reducing {name}: {stats_path}", flush=True)
        manifest["cases"].append(
            build_case(name, stats_path, sam_run, args.output_dir)
        )
    if not manifest["cases"]:
        raise SystemExit("No PDF-9 cases could be reduced")
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {args.output_dir / 'manifest.json'}")


if __name__ == "__main__":
    main()
