#!/usr/bin/env python3
"""Split one stats NetCDF file into legacy-style stats files.

Example:
  python utilities/split_stats_to_legacy.py \
    --input output/rico_silhs_stats.nc

Default output location is ./output, with legacy-style filenames, e.g.:
  output/rico_silhs_zt.nc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    import numpy as np
    from netCDF4 import Dataset
except Exception as exc:  # pragma: no cover - runtime dependency check
    print(f"ERROR: missing dependency: {exc}", file=sys.stderr)
    print("Install with: pip install netCDF4 numpy", file=sys.stderr)
    raise SystemExit(2)


COORD_AND_META_VARS = {
    "time",
    "col",
    "zt",
    "lh_zt",
    "zm",
    "rad_zt",
    "rad_zm",
    "lh_sample_number",
    "param",
    "param_name",
    "clubb_params",
}

OUTPUT_SUFFIXES = (
    "lh_zt",
    "sfc",
    "zm",
    "lh_sfc",
    "zt",
    "rad_zt",
    "rad_zm",
    "nl_lh_sample_points_2D",
    "u_lh_sample_points_2D",
)

LEGACY_LH_ZT_NAMES = {
    # Legacy files store these in *_lh_zt.nc even though they do not use an lh_ prefix.
    "AKm",
    "AKstd",
    "AKstd_cld",
    "AKm_rcm",
    "AKm_rcc",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split one *_stats.nc file into legacy-style stats files."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input stats file (e.g. output/rico_silhs_stats.nc)",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Output directory (default: output)",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Output file prefix (default: input stem with '_stats' removed)",
    )
    parser.add_argument(
        "--col-index",
        type=int,
        default=0,
        help="Column index to extract from the stats col dimension (default: 0)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files",
    )
    return parser.parse_args()


def derive_prefix(input_path: Path, explicit_prefix: str) -> str:
    if explicit_prefix:
        return explicit_prefix
    stem = input_path.stem
    if stem.endswith("_stats"):
        return stem[: -len("_stats")]
    return stem


def is_lh_like(name: str) -> bool:
    # Keep this heuristic simple; fields can be remapped later if desired.
    if name in LEGACY_LH_ZT_NAMES:
        return True
    if name.startswith("lh_"):
        return True
    if name.startswith("silhs_"):
        return True
    if name == "k_lh_start":
        return True
    return False


def classify_var(name: str, dims: Sequence[str]) -> str | None:
    dim_set = set(dims)

    if name.startswith("lh_nl_"):
        if {"time", "col", "lh_sample_number"}.issubset(dim_set) and (
            "zt" in dim_set or "lh_zt" in dim_set
        ):
            return "nl_lh_sample_points_2D"
        return None

    if name.startswith("lh_u_"):
        if {"time", "col", "lh_sample_number"}.issubset(dim_set) and (
            "zt" in dim_set or "lh_zt" in dim_set
        ):
            return "u_lh_sample_points_2D"
        return None

    if "time" not in dim_set or "col" not in dim_set:
        return None

    if "rad_zm" in dim_set:
        return "rad_zm"

    if "rad_zt" in dim_set:
        return "rad_zt"

    if "zm" in dim_set:
        return "zm"

    if "lh_zt" in dim_set:
        return "lh_zt"

    if "zt" in dim_set:
        return "lh_zt" if is_lh_like(name) else "zt"

    if len(dims) == 2:
        return "lh_sfc" if is_lh_like(name) else "sfc"

    return None


def copy_global_attrs(src: Dataset, dst: Dataset) -> None:
    for attr in src.ncattrs():
        dst.setncattr(attr, src.getncattr(attr))

    if "Conventions" not in dst.ncattrs():
        dst.setncattr("Conventions", "COARDS")
    if "model" not in dst.ncattrs():
        dst.setncattr("model", "CLUBB")


def init_legacy_file(
    out_path: Path,
    altitude_values: np.ndarray,
    time_values: np.ndarray,
    time_units: str,
    src: Dataset,
    overwrite: bool,
) -> Dataset:
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"{out_path} already exists (use --overwrite to replace)")

    ds = Dataset(out_path, "w", format="NETCDF4_CLASSIC")
    copy_global_attrs(src, ds)

    ds.createDimension("longitude", 1)
    ds.createDimension("latitude", 1)
    ds.createDimension("altitude", altitude_values.size)
    ds.createDimension("time", None)

    lon = ds.createVariable("longitude", "f8", ("longitude",))
    lon[:] = np.array([0.0], dtype=np.float64)
    lon.units = "degrees_E"

    lat = ds.createVariable("latitude", "f8", ("latitude",))
    lat[:] = np.array([0.0], dtype=np.float64)
    lat.units = "degrees_N"

    alt = ds.createVariable("altitude", "f8", ("altitude",))
    alt[:] = altitude_values
    alt.units = "meters"

    time_var = ds.createVariable("time", "f8", ("time",))
    time_var[:] = time_values
    time_var.units = time_units

    return ds


def init_lh_2d_file(
    out_path: Path,
    altitude_values: np.ndarray,
    time_values: np.ndarray,
    time_units: str,
    sample_values: np.ndarray,
    src: Dataset,
    overwrite: bool,
) -> Dataset:
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"{out_path} already exists (use --overwrite to replace)")

    ds = Dataset(out_path, "w", format="NETCDF4_CLASSIC")
    copy_global_attrs(src, ds)

    ds.createDimension("lh_sample_number", sample_values.size)
    ds.createDimension("longitude", 1)
    ds.createDimension("latitude", 1)
    ds.createDimension("altitude", altitude_values.size)
    ds.createDimension("time", None)

    sample = ds.createVariable("lh_sample_number", "f8", ("lh_sample_number",))
    sample[:] = sample_values
    sample.description = "SILHS sample (i.e. subcolumn) index"
    sample.units = "number"

    lon = ds.createVariable("longitude", "f8", ("longitude",))
    lon[:] = np.array([0.0], dtype=np.float64)
    lon.units = "degrees_E"

    lat = ds.createVariable("latitude", "f8", ("latitude",))
    lat[:] = np.array([0.0], dtype=np.float64)
    lat.units = "degrees_N"

    alt = ds.createVariable("altitude", "f8", ("altitude",))
    alt[:] = altitude_values
    alt.units = "meters"

    time_var = ds.createVariable("time", "f8", ("time",))
    time_var[:] = time_values
    time_var.units = time_units

    return ds


def move_col_to_last(arr: np.ndarray, dims: Sequence[str]) -> np.ndarray:
    col_idx = dims.index("col")
    if col_idx == arr.ndim - 1:
        return arr
    return np.moveaxis(arr, col_idx, -1)


def move_time_to_first(arr: np.ndarray, dims: Sequence[str]) -> np.ndarray:
    time_idx = dims.index("time")
    if time_idx == 0:
        return arr
    return np.moveaxis(arr, time_idx, 0)


def extract_legacy_payload(
    values: np.ndarray,
    dims: Sequence[str],
    col_index: int,
) -> np.ndarray:
    arr = np.asarray(values)
    arr = move_time_to_first(arr, dims)

    # Keep dims in sync after move_time_to_first.
    dims_after_time = list(dims)
    time_idx = dims_after_time.index("time")
    dims_after_time.pop(time_idx)
    dims_after_time.insert(0, "time")

    arr = move_col_to_last(arr, dims_after_time)

    dims_after_col = list(dims_after_time)
    col_idx = dims_after_col.index("col")
    dims_after_col.pop(col_idx)
    dims_after_col.append("col")

    if col_index < 0 or col_index >= arr.shape[-1]:
        raise IndexError(f"col-index {col_index} out of range for shape {arr.shape}")

    arr = arr[..., col_index]

    if arr.ndim == 1:
        # (time,) -> (time, altitude=1, latitude=1, longitude=1)
        return arr[:, np.newaxis, np.newaxis, np.newaxis]

    if arr.ndim == 2:
        # (time, altitude) -> (time, altitude, latitude=1, longitude=1)
        return arr[:, :, np.newaxis, np.newaxis]

    raise ValueError(f"Unsupported variable rank after selecting col: {arr.ndim}")


def copy_variable(
    src_var,
    dst: Dataset,
    col_index: int,
) -> None:
    # Do not copy _FillValue from stats: many fields use 0.0 there, and
    # legacy files do not mask zeros on read.
    dst_var = dst.createVariable(src_var.name, src_var.dtype, ("time", "altitude", "latitude", "longitude"))

    for attr in src_var.ncattrs():
        if attr == "_FillValue":
            continue
        dst_var.setncattr(attr, src_var.getncattr(attr))

    payload = extract_legacy_payload(src_var[:], src_var.dimensions, col_index)
    dst_var[:] = payload


def extract_lh_2d_payload(
    values: np.ndarray,
    dims: Sequence[str],
    col_index: int,
) -> np.ndarray:
    arr = np.asarray(values)
    arr = move_time_to_first(arr, dims)

    dims_after_time = list(dims)
    time_idx = dims_after_time.index("time")
    dims_after_time.pop(time_idx)
    dims_after_time.insert(0, "time")

    arr = move_col_to_last(arr, dims_after_time)

    dims_after_col = list(dims_after_time)
    col_idx = dims_after_col.index("col")
    dims_after_col.pop(col_idx)
    dims_after_col.append("col")

    if col_index < 0 or col_index >= arr.shape[-1]:
        raise IndexError(f"col-index {col_index} out of range for shape {arr.shape}")

    arr = arr[..., col_index]

    # Remaining dims should be exactly: time, {zt|lh_zt}, lh_sample_number.
    # This supports both LH source layouts:
    #   old: (lh_sample_number, col, lh_zt, time)
    #   new: (col, lh_sample_number, lh_zt, time)
    dims_after_select = [d for d in dims_after_col if d != "col"]

    if arr.ndim != len(dims_after_select) or arr.ndim != 3:
        raise ValueError(f"Unsupported LH 2D rank after selecting col: {arr.ndim}")

    time_axis = dims_after_select.index("time")
    if "lh_zt" in dims_after_select:
        alt_axis = dims_after_select.index("lh_zt")
    else:
        alt_axis = dims_after_select.index("zt")
    samp_axis = dims_after_select.index("lh_sample_number")
    arr = np.moveaxis(arr, (time_axis, alt_axis, samp_axis), (0, 1, 2))
    return arr[:, :, np.newaxis, np.newaxis, :]


def copy_variable_lh_2d(
    src_var,
    dst: Dataset,
    col_index: int,
    dst_name: str,
) -> None:
    dst_var = dst.createVariable(
        dst_name,
        src_var.dtype,
        ("time", "altitude", "latitude", "longitude", "lh_sample_number"),
    )

    for attr in src_var.ncattrs():
        if attr == "_FillValue":
            continue
        dst_var.setncattr(attr, src_var.getncattr(attr))

    payload = extract_lh_2d_payload(src_var[:], src_var.dimensions, col_index)
    dst_var[:] = payload


def collect_output_paths(output_dir: Path, prefix: str) -> Dict[str, Path]:
    return {name: output_dir / f"{prefix}_{name}.nc" for name in OUTPUT_SUFFIXES}


def main() -> int:
    args = parse_args()
    in_path = Path(args.input).resolve()
    if not in_path.exists():
        print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
        return 2

    output_dir = Path(args.output_dir).resolve()
    prefix = derive_prefix(in_path, args.prefix)
    out_paths = collect_output_paths(output_dir, prefix)

    output_dir.mkdir(parents=True, exist_ok=True)

    with Dataset(in_path) as src:
        # Read raw numeric values from stats. This avoids netCDF4 masking
        # all zeros when _FillValue=0.0 is present.
        src.set_auto_mask(False)

        if "time" not in src.variables:
            print("ERROR: input missing required coordinate variable 'time'", file=sys.stderr)
            return 2
        if "col" not in src.dimensions:
            print("ERROR: input missing required dimension 'col'", file=sys.stderr)
            return 2

        time_vals = np.asarray(src.variables["time"][:], dtype=np.float64)
        time_units = getattr(src.variables["time"], "units", "seconds")

        zt_vals = np.asarray(src.variables["zt"][:], dtype=np.float64) if "zt" in src.variables else np.array([], dtype=np.float64)
        lh_zt_vals = (
            np.asarray(src.variables["lh_zt"][:], dtype=np.float64)
            if "lh_zt" in src.variables
            else np.array([], dtype=np.float64)
        )
        zm_vals = np.asarray(src.variables["zm"][:], dtype=np.float64) if "zm" in src.variables else np.array([], dtype=np.float64)
        rad_zt_vals = np.asarray(src.variables["rad_zt"][:], dtype=np.float64) if "rad_zt" in src.variables else np.array([], dtype=np.float64)
        rad_zm_vals = np.asarray(src.variables["rad_zm"][:], dtype=np.float64) if "rad_zm" in src.variables else np.array([], dtype=np.float64)
        if "lh_sample_number" in src.variables:
            sample_vals = np.asarray(src.variables["lh_sample_number"][:], dtype=np.float64)
        elif "lh_sample_number" in src.dimensions:
            nsamp = len(src.dimensions["lh_sample_number"])
            sample_vals = np.asarray(np.arange(1, nsamp + 1), dtype=np.float64)
        else:
            sample_vals = np.array([], dtype=np.float64)
        sfc_vals = np.array([0.0], dtype=np.float64)

        altitudes = {
            "zt": zt_vals,
            "lh_zt": lh_zt_vals if lh_zt_vals.size > 0 else zt_vals,
            "zm": zm_vals,
            "sfc": sfc_vals,
            "lh_sfc": sfc_vals,
            "rad_zt": rad_zt_vals,
            "rad_zm": rad_zm_vals,
            "nl_lh_sample_points_2D": lh_zt_vals if lh_zt_vals.size > 0 else zt_vals,
            "u_lh_sample_points_2D": lh_zt_vals if lh_zt_vals.size > 0 else zt_vals,
        }

        dst_map: Dict[str, Dataset] = {}
        skipped_groups: List[str] = []
        for group, path in out_paths.items():
            # Do not emit profile files if their vertical grid is unavailable
            # in the input (e.g., rad_zt/rad_zm when stats lacked rad dims).
            if group not in ("sfc", "lh_sfc") and altitudes[group].size == 0:
                skipped_groups.append(group)
                continue
            if group in ("nl_lh_sample_points_2D", "u_lh_sample_points_2D"):
                if sample_vals.size == 0:
                    skipped_groups.append(group)
                    continue
                dst_map[group] = init_lh_2d_file(
                    out_path=path,
                    altitude_values=altitudes[group],
                    time_values=time_vals,
                    time_units=time_units,
                    sample_values=sample_vals,
                    src=src,
                    overwrite=args.overwrite,
                )
            else:
                dst_map[group] = init_legacy_file(
                    out_path=path,
                    altitude_values=altitudes[group],
                    time_values=time_vals,
                    time_units=time_units,
                    src=src,
                    overwrite=args.overwrite,
                )

        counts = {k: 0 for k in dst_map}
        skipped: List[Tuple[str, Tuple[str, ...]]] = []

        try:
            for name, var in src.variables.items():
                if name in COORD_AND_META_VARS:
                    continue

                group = classify_var(name, var.dimensions)
                if group is None:
                    skipped.append((name, tuple(var.dimensions)))
                    continue

                if group not in dst_map:
                    skipped.append((name, tuple(var.dimensions)))
                    continue

                if group == "nl_lh_sample_points_2D":
                    copy_variable_lh_2d(var, dst_map[group], args.col_index, name[len("lh_nl_"):])
                elif group == "u_lh_sample_points_2D":
                    copy_variable_lh_2d(var, dst_map[group], args.col_index, name[len("lh_u_"):])
                else:
                    copy_variable(var, dst_map[group], args.col_index)
                counts[group] += 1
        finally:
            for ds in dst_map.values():
                ds.close()

    print(f"Wrote files with prefix '{prefix}' in {output_dir}")
    for group in OUTPUT_SUFFIXES:
        if group in dst_map:
            print(f"  - {out_paths[group].name}: {counts[group]} vars")
    if skipped_groups:
        print("Skipped output files due to missing vertical grid in input:")
        for group in skipped_groups:
            print(f"  - {out_paths[group].name}")
    if skipped:
        print(f"Skipped {len(skipped)} vars (unsupported dims):")
        for name, dims in skipped[:20]:
            print(f"  - {name}: {dims}")
        if len(skipped) > 20:
            print(f"  ... {len(skipped) - 20} more")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
