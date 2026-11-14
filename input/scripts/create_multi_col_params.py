#!/usr/bin/env python3
"""Build multicolumn parameter arrays for runs and explicit JSON exports.

Usage:
    multi_col_settings() transforms resolved settings in memory. The CLI reads a
    catalog or sparse parameter file and writes the same result to -output_file.
    is_single_col_list() checks a file's dimensions; parse_clubb_params() reads
    its scalar model values. duplicate_params() copies values, while
    duplicate_and_tweak() varies selected parameters and optionally mirrors them.
    hypergrid() samples relative ranges. parse_hypergrid_range_spec() reads
    explicit ranges for custom_hypergrid(). clamp_batch_size() caps batch width.

Understanding:
    Generated lists contain one value per column in meshgrid order. SILHS
    settings are copied unchanged. Input must be single-column; the default
    count mode varies C11 over half to twice its initial value. Exact duplication
    is used for tuning batches. Batch width is independent of the column count.
"""
import argparse
import json
import sys
import os
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

import numpy as np

from input.scripts.parameters_json_to_namelist import parameter_settings, validate_parameter_settings
from input.scripts.namelist import load_json


def _parse_float(text):
    """Read either ordinary or Fortran-style exponent notation in parameter grids."""
    return float(text.replace("D", "E").replace("d", "e"))


def _parameter_file_settings(path):
    """Resolve a catalog's defaults or a custom parameter override file."""
    data = load_json(path)
    return parameter_settings(json_file=path) if isinstance(data, dict) and "presets" in data else parameter_settings(path)


def is_single_col_list(path: str) -> bool:
    """Return whether a parameter JSON snapshot defines a single column."""
    return "multicol" not in _parameter_file_settings(path)


def parse_clubb_params(file_name):
    """Return the model parameter values from the unified parameter JSON."""
    return {name: float(value) for name, value in _parameter_file_settings(file_name)["model"].items()}


# Copy the parameters stored in parsed_params ngrdcol times
# For exampe:
#   C1 = 1.0
# becomes
#   C1 = 1.0, 1.0, 1.0, ...   <- ngrdcol entries
def duplicate_params( ngrdcol, parsed_params ):

    clubb_params = {}

    for var_name, value in parsed_params.items():
        clubb_params[var_name] = np.tile(value, ngrdcol)

    return clubb_params


# Duplicate the parameters, then change them slightly based
# on a predefined list and tweaking behavior
def duplicate_and_tweak( ngrdcol, parsed_params, tweak_list, mirror=False ):

    clubb_params = duplicate_params( ngrdcol, parsed_params )

    for param in tweak_list:
        # Tweak parameters in the list to have ngrdcol sample points
        # in the range initial_param_val*[1/2, 2]
        initial_param_val   = parsed_params[param]
        clubb_params[param] = initial_param_val * np.linspace( 0.5, 2, ngrdcol)

        # If mirror is set we mirror the first half of the parameter set to the last half
        # This is mainly for testing purposes
        if mirror:
            first_half = clubb_params[param][:ngrdcol//2]
            clubb_params[param][(ngrdcol+1)//2:] = first_half[::-1]

    return clubb_params


# Create a hypergrid from a list a variables to tweak
# For exmaple:
#   -tweak_list C7,C8  and  -ngrdcol 10
# will result in a grid array of points with each parameter 
# axis being 10 possible values
def hypergrid(samples_per_param, parsed_params, tweak_list):
    ngrdcol = samples_per_param ** len(tweak_list)

    # Create points for each parameter in tweak_list
    param_grids = []
    for param in tweak_list:
        param_points = parsed_params[param] * np.linspace(0.5, 2, samples_per_param)
        param_grids.append(param_points)

    # Generate a meshgrid for all parameters
    mesh = np.meshgrid(*param_grids, indexing='ij')

    # Flatten each parameter's grid and store in clubb_params
    clubb_params = duplicate_params(ngrdcol, parsed_params)
    for i, param in enumerate(tweak_list):
        clubb_params[param][:] = mesh[i].ravel()

    return clubb_params, ngrdcol


def parse_hypergrid_range_spec(spec_text):
    """
    Parse a custom hypergrid specification of the form:
      PARAM1/MIN:MAX/NPOINTS,PARAM2/MIN:MAX/NPOINTS,...

    One coordinate may deliberately target an equal-value parameter group:
      C6rt=C6thl/0:4/5
    The group still consumes one hypergrid axis and assigns the same sampled
    value to every target in each column.
    """
    specs = []
    for raw_item in (spec_text or "").split(","):
        item = raw_item.strip()
        if not item:
            continue
        parts = item.split("/")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid -multicol entry '{item}'. Expected PARAM/MIN:MAX/NPOINTS."
            )
        coordinate_name = parts[0].strip()
        targets = tuple(name.strip() for name in coordinate_name.split("="))
        if not coordinate_name or not targets or any(not name for name in targets):
            raise ValueError(f"Invalid -multicol entry '{item}': missing parameter name.")
        try:
            min_text, max_text = parts[1].split(":", 1)
        except ValueError as exc:
            raise ValueError(
                f"Invalid -multicol entry '{item}'. Range must look like MIN:MAX."
            ) from exc
        min_value = _parse_float(min_text.strip())
        max_value = _parse_float(max_text.strip())
        npoints = int(parts[2].strip())
        if npoints < 1:
            raise ValueError(
                f"Invalid -multicol entry '{item}'. NPOINTS must be >= 1."
            )
        specs.append(
            {
                "name": "=".join(targets),
                "targets": targets,
                "min": min_value,
                "max": max_value,
                "npoints": npoints,
            }
        )
    if not specs:
        raise ValueError("The -multicol specification did not contain any parameter ranges.")
    return specs


def custom_hypergrid(parsed_params, range_specs):
    """Create a hypergrid using explicit per-parameter ranges and point counts."""
    ngrdcol = 1
    param_grids = []
    coordinates = []
    claimed_targets = set()
    for spec in range_specs:
        targets = tuple(spec.get("targets") or (spec["name"],))
        duplicate_targets = claimed_targets.intersection(targets)
        if duplicate_targets:
            names = ", ".join(sorted(duplicate_targets))
            raise ValueError(f"Parameter target(s) requested by more than one hypergrid coordinate: {names}.")
        for param_name in targets:
            if param_name not in parsed_params:
                raise KeyError(f"Parameter '{param_name}' was not found in the input params file.")
        claimed_targets.update(targets)
        coordinates.append(targets)
        npoints = spec["npoints"]
        ngrdcol *= npoints
        param_grids.append(np.linspace(spec["min"], spec["max"], npoints))

    mesh = np.meshgrid(*param_grids, indexing='ij')
    clubb_params = duplicate_params(ngrdcol, parsed_params)
    for i, targets in enumerate(coordinates):
        values = mesh[i].ravel()
        for param_name in targets:
            clubb_params[param_name][:] = values

    return clubb_params, ngrdcol


def clamp_batch_size(batch_size, ngrdcol):
    """Keep runtime batch size within the generated column count."""
    if batch_size is None:
        return None
    if batch_size <= ngrdcol:
        return batch_size
    print(
        f" - requested batch_size = {batch_size} exceeds ngrdcol = {ngrdcol}; "
        f"using batch_size = {ngrdcol}"
    )
    return ngrdcol


def multi_col_settings(settings, multicol, *, batch_size=None, mode="dup_tweak",
                       tweak_list=None, mirror=False):
    """Copy resolved single-column settings into arrays, retaining SILHS values.

    multicol is a count or a PARAM/MIN:MAX/NPOINTS grid specification. Count
    modes duplicate, duplicate and tweak, or form a relative hypergrid. Return
    plain JSON values with ngrdcol and an optional clipped batch_size.
    """
    import copy

    if "multicol" in settings:
        raise ValueError("Column generation requires single-column parameters")
    if batch_size is not None and batch_size < 1:
        raise ValueError("-batch_size must be >= 1")
    parsed = {name: float(value) for name, value in settings["model"].items()}
    tweaks = ["C11"] if tweak_list is None else tweak_list
    try:
        count = int(multicol)
    except ValueError:
        values, count = custom_hypergrid(parsed, parse_hypergrid_range_spec(multicol))
    else:
        if count < 1:
            raise ValueError("-multicol column count must be >= 1")
        if mode == "duplicate":
            values = duplicate_params(count, parsed)
        elif mode == "dup_tweak":
            values = duplicate_and_tweak(count, parsed, tweaks, mirror)
        elif mode == "hypergrid":
            values, count = hypergrid(count, parsed, tweaks)
        else:
            raise ValueError(f"Unknown column generation mode: {mode}")
    result = copy.deepcopy(settings)
    result["model"] = {name: value.tolist() for name, value in values.items()}
    result["multicol"] = {"ngrdcol": count}
    if batch_size is not None:
        result["multicol"]["batch_size"] = clamp_batch_size(batch_size, count)
    return validate_parameter_settings(result)


if __name__ == "__main__":

    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process CLUBB parameters.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")

    parser.add_argument("-multicol", required=True, metavar="NUM|SPEC",
                        help="Column count or PARAM/MIN:MAX/NPOINTS parameter grid.")

    tunable_parameters_default_path = os.path.join(script_dir, "../clubb_params.jsonc")
    parser.add_argument( '-params_file', dest='param_file', type=str, help="Path to the CLUBB parameters file",
                         default = str(tunable_parameters_default_path) , metavar='FILE')

    parser.add_argument( '-output_file', dest='out_file', type=str, help="Output parameter JSON snapshot",
                         default = "clubb_params_multi_col.json" , metavar='FILE')

    parser.add_argument( "-mode", type=str, choices=("duplicate", "dup_tweak", "hypergrid"),
                         help="Count mode (default: dup_tweak); hypergrid interprets NUM as points per axis." )

    parser.add_argument(
        "-batch_size",
        type=int,
        help=(
            "Runtime batch size written to &multicol_def. "
            "Values larger than ngrdcol are clipped to ngrdcol."
        ),
    )

    parser.add_argument( "-mirror", type=str, help="mirror param lists",
                         default = "false" )

    parser.add_argument(
        "-tweak_list",
        type=lambda s: s.split(','),  # Split the input string on commas
        help="Comma-separated list of tweaks (e.g., C1,C2,C3)"
    )


    # Parse the arguments
    args = parser.parse_args()

    try:
        int(args.multicol)
    except ValueError:
        if args.mode is not None:
            parser.error("Do not specify -mode when using a parameter-grid specification.")

    try:
        settings = multi_col_settings(
            _parameter_file_settings(args.param_file), args.multicol,
            batch_size=args.batch_size, mode=args.mode or "dup_tweak",
            tweak_list=args.tweak_list, mirror=args.mirror == "true",
        )
    except (KeyError, ValueError) as exc:
        parser.error(str(exc))
    from input.scripts.namelist import atomic_text
    atomic_text(args.out_file, json.dumps(settings, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {settings['multicol']['ngrdcol']} columns to '{args.out_file}'")
