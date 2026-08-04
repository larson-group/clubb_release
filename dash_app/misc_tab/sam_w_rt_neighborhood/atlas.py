#!/usr/bin/env python3
"""Generate static sprite atlases for the raw-SAM w-r_t Misc explorer."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
import time

import netCDF4
import numpy as np
from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[3]

from dash_app.shared.bivariate_heatmap import signed_bivariate_reference_rgb
from utilities.sam_3d_reference import DEFAULT_SAM_RUN


ATLAS_VERSION = 1
STEP_PATTERN = re.compile(r"_(\d+)_micro\.nc$")


def _snapshot_files(run_dir: Path) -> list[tuple[int, Path]]:
    pairs = []
    for path in (run_dir / "OUT_3D").glob("*_micro.nc"):
        match = STEP_PATTERN.search(path.name)
        if match:
            pairs.append((int(match.group(1)), path))
    pairs.sort()
    if not pairs:
        raise FileNotFoundError(f"No SAM *_micro.nc files found under {run_dir / 'OUT_3D'}")
    return pairs


def _read_fields(path: Path, spatial_stride: int = 1):
    with netCDF4.Dataset(path) as dataset:
        slicing = (0, slice(None), slice(None, None, spatial_stride), slice(None, None, spatial_stride))
        w = np.asarray(dataset.variables["W"][slicing], dtype=np.float64)
        rt = np.asarray(dataset.variables["RT"][slicing], dtype=np.float64) * 1000.0
        rc = np.asarray(dataset.variables["RC"][slicing], dtype=np.float64)
        heights = np.asarray(dataset.variables["z"][:], dtype=np.float64)
    return (
        w.reshape(w.shape[0], -1),
        rt.reshape(rt.shape[0], -1),
        rc.reshape(rc.shape[0], -1),
        heights,
    )


def _sample_axis_ranges(sample_fields, height_count: int):
    x_ranges = np.empty((height_count, 2), dtype=np.float64)
    y_ranges = np.empty((height_count, 2), dtype=np.float64)
    for level in range(height_count):
        w = np.concatenate([fields[0][level] for fields in sample_fields])
        rt = np.concatenate([fields[1][level] for fields in sample_fields])
        rc = np.concatenate([fields[2][level] for fields in sample_fields])
        x_low, x_high = np.percentile(rt, (0.05, 99.95))
        y_low, y_high = np.percentile(w, (0.05, 99.95))
        cloudy = rc > 0.0
        if np.any(cloudy):
            cloud_x_low, cloud_x_high = np.percentile(rt[cloudy], (0.5, 99.5))
            cloud_y_low, cloud_y_high = np.percentile(w[cloudy], (0.5, 99.5))
            x_low = min(x_low, cloud_x_low)
            x_high = max(x_high, cloud_x_high)
            y_low = min(y_low, cloud_y_low)
            y_high = max(y_high, cloud_y_high)
        x_pad = max(0.04 * (x_high - x_low), 1.0e-3)
        y_pad = max(0.04 * (y_high - y_low), 1.0e-3)
        x_ranges[level] = (x_low - x_pad, x_high + x_pad)
        y_ranges[level] = (y_low - y_pad, y_high + y_pad)
    return x_ranges, y_ranges


def _binned_fields(w, rt, rc, x_ranges, y_ranges, bins: int):
    height_count, sample_count = w.shape
    x_scale = bins / np.maximum(x_ranges[:, 1] - x_ranges[:, 0], 1.0e-12)
    y_scale = bins / np.maximum(y_ranges[:, 1] - y_ranges[:, 0], 1.0e-12)
    x_bin = np.floor((rt - x_ranges[:, 0, None]) * x_scale[:, None]).astype(np.int64)
    y_bin = np.floor((w - y_ranges[:, 0, None]) * y_scale[:, None]).astype(np.int64)
    valid = (x_bin >= 0) & (x_bin < bins) & (y_bin >= 0) & (y_bin < bins)
    level_index = np.broadcast_to(np.arange(height_count)[:, None], w.shape)
    flat_index = level_index * bins * bins + y_bin * bins + x_bin
    probability = np.bincount(
        flat_index[valid], minlength=height_count * bins * bins
    ).reshape(height_count, bins, bins)
    probability = probability.astype(np.float64) / float(sample_count)

    centered_w = w - np.mean(w, axis=1, keepdims=True)
    centered_rc = rc - np.mean(rc, axis=1, keepdims=True)
    signed_local = centered_w * centered_rc * 1000.0
    signed = np.bincount(
        flat_index[valid],
        weights=signed_local[valid],
        minlength=height_count * bins * bins,
    ).reshape(height_count, bins, bins)
    signed /= float(sample_count)
    return probability, signed


def _panel_scale(values: np.ndarray) -> np.ndarray:
    scales = np.zeros(values.shape[0], dtype=np.float64)
    for level, field in enumerate(values):
        positive = np.asarray(field, dtype=np.float64)
        positive = positive[np.isfinite(positive) & (positive > 0.0)]
        if positive.size:
            scales[level] = float(np.percentile(positive, 99.5))
    return scales


def _sample_color_scales(sample_fields, x_ranges, y_ranges, bins: int):
    probability_scales = []
    transport_scales = []
    for w, rt, rc, _heights in sample_fields:
        probability, signed = _binned_fields(w, rt, rc, x_ranges, y_ranges, bins)
        probability_scales.append(_panel_scale(probability))
        transport_scales.append(_panel_scale(np.abs(signed)))
    probability_upper = np.percentile(np.stack(probability_scales), 95.0, axis=0)
    transport_upper = np.percentile(np.stack(transport_scales), 95.0, axis=0)
    probability_upper = np.maximum(probability_upper, np.finfo(np.float64).tiny)
    transport_upper = np.maximum(transport_upper, np.finfo(np.float64).tiny)
    return probability_upper, transport_upper


def _normalize(field, upper):
    field = np.nan_to_num(field, nan=0.0, posinf=0.0, neginf=0.0)
    upper = np.maximum(np.asarray(upper, dtype=np.float64), np.finfo(np.float64).tiny)
    knee = np.maximum(0.025 * upper, np.finfo(np.float64).tiny)
    normalized = np.log1p(np.maximum(field, 0.0) / knee[:, None, None])
    normalized /= np.log1p(upper / knee)[:, None, None]
    return np.clip(normalized, 0.0, 1.0)


def _color_codes(probability, signed, probability_upper, transport_upper, levels: int):
    normalized_probability = _normalize(probability, probability_upper)
    normalized_magnitude = _normalize(np.abs(signed), transport_upper)
    normalized_signed = np.sign(signed) * normalized_magnitude
    probability_bin = np.rint(normalized_probability * (levels - 1)).astype(np.uint8)
    signed_bin = np.rint(normalized_signed * (levels - 1)).astype(np.int16) + levels - 1
    return (probability_bin.astype(np.int16) * (2 * levels - 1) + signed_bin).astype(np.uint8)


def _palette(levels: int) -> list[int]:
    colors = np.zeros((256, 3), dtype=np.uint8)
    signed_levels = 2 * levels - 1
    for probability_bin in range(levels):
        for signed_bin in range(signed_levels):
            code = probability_bin * signed_levels + signed_bin
            probability = probability_bin / (levels - 1)
            signed_transport = (signed_bin - (levels - 1)) / (levels - 1)
            colors[code] = signed_bivariate_reference_rgb(probability, signed_transport)
    return colors.reshape(-1).tolist()


def _write_manifest(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def generate_atlas(
    run_dir: Path,
    output_dir: Path,
    *,
    bins: int = 48,
    time_block: int = 20,
    height_block: int = 10,
    levels: int = 8,
    axis_time_stride: int = 10,
    axis_spatial_stride: int = 4,
    force: bool = False,
) -> Path:
    """Generate sprite atlases and return the completed manifest path."""

    run_dir = Path(run_dir).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    pairs = _snapshot_files(run_dir)
    if levels * (2 * levels - 1) > 256:
        raise ValueError("Palette atlases require levels*(2*levels-1) <= 256.")
    if bins < 8 or time_block < 5 or height_block < 5:
        raise ValueError("bins must be >=8 and atlas blocks must each contain at least 5 cells.")
    if output_dir.exists() and not force:
        manifest = output_dir / "manifest.json"
        if manifest.exists():
            return manifest
        raise FileExistsError(f"Incomplete output exists at {output_dir}; use --force.")

    sample_indices = list(range(0, len(pairs), max(int(axis_time_stride), 1)))
    if sample_indices[-1] != len(pairs) - 1:
        sample_indices.append(len(pairs) - 1)
    print(f"Sampling {len(sample_indices)} of {len(pairs)} times for stable axes and colors…", flush=True)
    sample_fields = [
        _read_fields(pairs[index][1], max(int(axis_spatial_stride), 1))
        for index in sample_indices
    ]
    heights = sample_fields[0][3]
    height_count = len(heights)
    if any(not np.array_equal(fields[3], heights) for fields in sample_fields[1:]):
        raise ValueError("SAM height coordinates change between sampled files.")
    x_ranges, y_ranges = _sample_axis_ranges(sample_fields, height_count)
    probability_upper, transport_upper = _sample_color_scales(
        sample_fields, x_ranges, y_ranges, bins
    )
    del sample_fields

    parent = output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}-", dir=parent))
    palette = _palette(levels)
    time_block_count = math.ceil(len(pairs) / time_block)
    height_block_count = math.ceil(height_count / height_block)
    started = time.perf_counter()
    try:
        for time_block_index in range(time_block_count):
            time_start = time_block_index * time_block
            atlas_arrays = [
                np.full(
                    (height_block * bins, time_block * bins),
                    levels - 1,
                    dtype=np.uint8,
                )
                for _ in range(height_block_count)
            ]
            for local_time in range(time_block):
                time_index = time_start + local_time
                if time_index >= len(pairs):
                    break
                w, rt, rc, current_heights = _read_fields(pairs[time_index][1])
                if not np.array_equal(current_heights, heights):
                    raise ValueError(f"Height coordinates changed in {pairs[time_index][1]}")
                probability, signed = _binned_fields(w, rt, rc, x_ranges, y_ranges, bins)
                codes = _color_codes(
                    probability,
                    signed,
                    probability_upper,
                    transport_upper,
                    levels,
                )
                codes = codes[:, ::-1, :]
                for height_index in range(height_count):
                    block_index = height_index // height_block
                    local_height = height_index % height_block
                    row = slice(local_height * bins, (local_height + 1) * bins)
                    col = slice(local_time * bins, (local_time + 1) * bins)
                    atlas_arrays[block_index][row, col] = codes[height_index]

            for height_block_index, atlas in enumerate(atlas_arrays):
                height_start = height_block_index * height_block
                image = Image.fromarray(atlas, mode="P")
                image.putpalette(palette)
                image.save(
                    temporary
                    / f"atlas_t{time_start:04d}_z{height_start:03d}.png",
                    format="PNG",
                    compress_level=6,
                )
            elapsed = time.perf_counter() - started
            print(
                f"Rendered time block {time_block_index + 1}/{time_block_count} "
                f"({min(time_start + time_block, len(pairs))}/{len(pairs)} times, {elapsed:.1f} s)",
                flush=True,
            )

        manifest_payload = {
            "version": ATLAS_VERSION,
            "projection": "w_rt",
            "coloring": "probability_plus_signed_wprcp",
            "source_run": str(run_dir),
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "time_seconds": [step for step, _path in pairs],
            "height_m": heights.tolist(),
            "tile_pixels": bins,
            "time_block": time_block,
            "height_block": height_block,
            "palette_levels": levels,
            "atlas_filename": "atlas_t{time_start:04d}_z{height_start:03d}.png",
            "x_range_rt_g_per_kg": x_ranges.tolist(),
            "y_range_w_m_per_s": y_ranges.tolist(),
            "probability_upper": probability_upper.tolist(),
            "signed_transport_abs_upper_m_per_s_g_per_kg": transport_upper.tolist(),
            "axis_sampling": {
                "time_stride": axis_time_stride,
                "horizontal_stride": axis_spatial_stride,
                "all_air_percentiles": [0.05, 99.95],
                "cloudy_air_percentiles": [0.5, 99.5],
            },
        }
        _write_manifest(temporary / "manifest.json", manifest_payload)
        if output_dir.exists():
            shutil.rmtree(output_dir)
        temporary.replace(output_dir)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(f"Wrote {output_dir / 'manifest.json'}", flush=True)
    return output_dir / "manifest.json"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sam-run", type=Path, default=DEFAULT_SAM_RUN)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--bins", type=int, default=48)
    parser.add_argument("--time-block", type=int, default=20)
    parser.add_argument("--height-block", type=int, default=10)
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("--axis-time-stride", type=int, default=10)
    parser.add_argument("--axis-spatial-stride", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = args.sam_run / "PDF_TILES" / "w_rt_signed"
    generate_atlas(
        args.sam_run,
        output_dir,
        bins=args.bins,
        time_block=args.time_block,
        height_block=args.height_block,
        levels=args.levels,
        axis_time_stride=args.axis_time_stride,
        axis_spatial_stride=args.axis_spatial_stride,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
