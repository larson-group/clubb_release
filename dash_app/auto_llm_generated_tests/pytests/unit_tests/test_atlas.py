"""Small end-to-end contract test for the static raw-SAM sprite generator."""

from __future__ import annotations

import json
import math
from pathlib import Path

import netCDF4
import numpy as np
from PIL import Image

from dash_app.misc_tab.sam_w_rt_neighborhood.atlas import generate_atlas


def _write_snapshot(path: Path, step: int, height_count: int = 7) -> None:
    y_count = x_count = 8
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("time", 1)
        dataset.createDimension("z", height_count)
        dataset.createDimension("y", y_count)
        dataset.createDimension("x", x_count)
        height = dataset.createVariable("z", "f8", ("z",))
        w_variable = dataset.createVariable("W", "f4", ("time", "z", "y", "x"))
        rt_variable = dataset.createVariable("RT", "f4", ("time", "z", "y", "x"))
        rc_variable = dataset.createVariable("RC", "f4", ("time", "z", "y", "x"))
        height[:] = 20.0 + 40.0 * np.arange(height_count)
        yy, xx = np.meshgrid(
            np.linspace(-1.0, 1.0, y_count),
            np.linspace(-1.0, 1.0, x_count),
            indexing="ij",
        )
        for level in range(height_count):
            phase = 0.04 * step + 0.15 * level
            w = 0.6 * xx + 0.25 * yy + phase
            rt = 0.012 + 0.0015 * xx + 0.0002 * level + 1.0e-7 * step
            rc = np.maximum(rt - (0.0128 + 0.00015 * level), 0.0)
            w_variable[0, level] = w
            rt_variable[0, level] = rt
            rc_variable[0, level] = rc


def test_generate_static_sprite_atlas(tmp_path):
    run_dir = tmp_path / "sam"
    for index in range(7):
        step = 60 * (index + 1)
        _write_snapshot(
            run_dir / "OUT_3D" / f"arm_{step:010d}_micro.nc",
            step,
        )

    output_dir = tmp_path / "atlas"
    manifest_path = generate_atlas(
        run_dir,
        output_dir,
        bins=8,
        time_block=5,
        height_block=5,
        levels=4,
        axis_time_stride=2,
        axis_spatial_stride=2,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["projection"] == "w_rt"
    assert len(manifest["time_seconds"]) == 7
    assert len(manifest["height_m"]) == 7

    images = sorted(output_dir.glob("atlas_*.png"))
    assert len(images) == math.ceil(7 / 5) ** 2
    with Image.open(images[0]) as image:
        assert image.mode == "P"
        assert image.size == (40, 40)
        assert np.unique(np.asarray(image)).size > 1

    # A complete atlas is deliberately reusable unless the caller asks to replace it.
    assert generate_atlas(run_dir, output_dir) == manifest_path
