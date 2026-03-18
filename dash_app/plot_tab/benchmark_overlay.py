import ast
import os
from collections import OrderedDict

import numpy as np
from netCDF4 import Dataset


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PYPLOTGEN_ROOT = os.path.join(REPO_ROOT, "postprocessing", "pyplotgen")

Z_DIM_NAMES = {"z", "zm", "zt", "altitude", "height", "lev"}
T_DIM_NAMES = {"t", "time"}
BENCHMARK_SOURCES = ("sam", "coamps")

_CACHE_MAX_ENTRIES = 256
_CATALOG_CACHE = OrderedDict()
_DATA_CACHE = OrderedDict()
_META_CACHE = OrderedDict()
_SOURCE_CACHE = OrderedDict()


def clear_benchmark_caches():
    _CATALOG_CACHE.clear()
    _DATA_CACHE.clear()
    _META_CACHE.clear()
    _SOURCE_CACHE.clear()


def available_benchmark_sources(case_data):
    benchmarks = (case_data or {}).get("benchmarks") or {}
    return list(benchmarks.get("available_sources") or [])


def sanitize_enabled_sources(case_data, enabled_sources):
    available = set(available_benchmark_sources(case_data))
    return [source for source in (enabled_sources or []) if source in available]


def _file_signature(path):
    try:
        stat = os.stat(path)
    except OSError:
        return (path, None, None)
    return (path, stat.st_mtime_ns, stat.st_size)


def _cached_lru(cache, max_entries, key, builder):
    cached = cache.get(key)
    if cached is not None:
        cache.move_to_end(key)
        return cached
    value = builder()
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > max_entries:
        cache.popitem(last=False)
    return value


def _freeze(value):
    if isinstance(value, np.ndarray):
        arr = np.asarray(value)
        arr.setflags(write=False)
        return arr
    if isinstance(value, dict):
        return {key: _freeze(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze(item) for item in value)
    return value


def _normalize_benchmark_paths(raw_value):
    if raw_value is None:
        return None
    if isinstance(raw_value, str):
        candidates = {"default": raw_value}
    elif isinstance(raw_value, dict):
        candidates = dict(raw_value)
    else:
        return None
    resolved = {}
    for key, value in candidates.items():
        if not value:
            continue
        path = os.path.abspath(os.path.expanduser(str(value)))
        if os.path.isfile(path):
            resolved[str(key)] = path
    return resolved or None


def build_case_benchmark_info(case_def):
    case_def = case_def or {}
    sam_paths = _normalize_benchmark_paths(case_def.get("sam_benchmark_file"))
    coamps_paths = _normalize_benchmark_paths(case_def.get("coamps_benchmark_file"))
    available = []
    if sam_paths:
        available.append("sam")
    if coamps_paths:
        available.append("coamps")
    return {
        "sam": sam_paths,
        "coamps": coamps_paths,
        "available_sources": available,
        "var_group_names": [str(name) for name in (case_def.get("var_groups") or []) if name],
    }


def _find_dim(dim_names, candidates):
    for dim_name in dim_names:
        if dim_name.lower() in candidates:
            return dim_name
    return None


def _time_units_factor(units):
    if not units:
        return 1.0
    unit = units.split()[0].lower()
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return 1.0
    if unit in {"min", "mins", "minute", "minutes"}:
        return 60.0
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return 3600.0
    return 1.0


def _time_values_seconds(values, units):
    arr = np.asarray(values, dtype=float)
    return arr * float(_time_units_factor(units))


def _slider_value_to_index(value, max_len):
    if max_len <= 0:
        return 0
    return max(0, min(int(value or 1) - 1, max_len - 1))


def _slider_range_to_indices(value, max_len):
    if max_len <= 0:
        return [0, 0]
    if not value or len(value) != 2:
        return [0, max_len - 1]
    start = _slider_value_to_index(value[0], max_len)
    end = _slider_value_to_index(value[1], max_len)
    if end < start:
        start, end = end, start
    return [start, end]


def _elapsed_seconds(values):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    return arr - float(arr[0])


def _profile_label(source_name):
    if source_name == "sam":
        return "SAM benchmark"
    if source_name == "coamps":
        return "COAMPS benchmark"
    return f"{source_name} benchmark"


def _benchmark_line_style(source_name):
    if source_name == "sam":
        return {"width": 4.0, "dash": "solid", "color": "#facc15"}
    if source_name == "coamps":
        return {"width": 4.0, "dash": "solid", "color": "#b3b3b3"}
    return {"width": 4.0, "dash": "solid"}


def benchmark_trace_style(source_name):
    return dict(_benchmark_line_style(source_name))


def _load_variable_definitions(group_name):
    def _build():
        path = os.path.join(PYPLOTGEN_ROOT, "config", f"{group_name}.py")
        if not os.path.isfile(path):
            return []
        with open(path, "r", encoding="utf-8") as handle:
            tree = ast.parse(handle.read(), filename=path)
        for node in tree.body:
            if not isinstance(node, ast.ClassDef) or node.name != group_name:
                continue
            for class_item in node.body:
                if not isinstance(class_item, ast.FunctionDef) or class_item.name != "__init__":
                    continue
                for stmt in class_item.body:
                    if not isinstance(stmt, ast.Assign):
                        continue
                    for target in stmt.targets:
                        if (
                            isinstance(target, ast.Attribute)
                            and isinstance(target.value, ast.Name)
                            and target.value.id == "self"
                            and target.attr == "variable_definitions"
                        ):
                            value = _eval_ast(stmt.value)
                            return list(value or []) if isinstance(value, list) else []
        return []

    return _cached_lru(_CATALOG_CACHE, _CACHE_MAX_ENTRIES, ("var_group_defs", group_name), _build)


def _all_variable_group_names():
    def _build():
        config_dir = os.path.join(PYPLOTGEN_ROOT, "config")
        names = []
        for filename in sorted(os.listdir(config_dir)):
            if not filename.startswith("VariableGroup") or not filename.endswith(".py"):
                continue
            if filename == "__init__.py":
                continue
            names.append(filename[:-3])
        return tuple(names)

    return _cached_lru(_CATALOG_CACHE, _CACHE_MAX_ENTRIES, ("all_var_groups",), _build)


def _eval_ast(node):
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.List):
        items = []
        for elt in node.elts:
            value = _eval_ast(elt)
            if value is not None:
                items.append(value)
        return items
    if isinstance(node, ast.Tuple):
        items = []
        for elt in node.elts:
            value = _eval_ast(elt)
            if value is not None:
                items.append(value)
        return tuple(items)
    if isinstance(node, ast.Dict):
        result = {}
        for key_node, value_node in zip(node.keys, node.values):
            key = _eval_ast(key_node)
            value = _eval_ast(value_node)
            if key is None or value is None:
                continue
            result[key] = value
        return result
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        value = _eval_ast(node.operand)
        if isinstance(value, (int, float)):
            return -value
        return None
    if isinstance(node, ast.BinOp):
        left = _eval_ast(node.left)
        right = _eval_ast(node.right)
        if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
            return None
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        return None
    if isinstance(node, ast.Name):
        return None
    if isinstance(node, ast.Attribute):
        if isinstance(node.value, ast.Name) and node.value.id == "self":
            return {"__call__": node.attr}
        return node.attr
    if isinstance(node, ast.Call):
        return None
    return None


def _build_benchmark_catalog(group_names):
        catalog = {}
        for group_name in group_names:
            for definition in _load_variable_definitions(group_name):
                if definition.get("type", "profile") not in {"profile", None}:
                    continue
                var_names = definition.get("var_names") or {}
                clubb_names = [name for name in (var_names.get("clubb") or []) if isinstance(name, str)]
                if not clubb_names:
                    continue
                for clubb_name in clubb_names:
                    entry = catalog.setdefault(
                        clubb_name,
                        {
                            "sam": {"entries": [], "conversion_factor": None},
                            "coamps": {"entries": [], "conversion_factor": None},
                        },
                    )
                    for source_name in BENCHMARK_SOURCES:
                        source_entry = entry[source_name]
                        for item in (var_names.get(source_name) or []):
                            mapped = None
                            if isinstance(item, str):
                                mapped = {"kind": "var", "name": item}
                            elif isinstance(item, dict) and "__call__" in item:
                                mapped = {"kind": "calc", "name": str(item["__call__"])}
                            if mapped is not None and mapped not in source_entry["entries"]:
                                source_entry["entries"].append(mapped)
                        conv_key = f"{source_name}_conv_factor"
                        if source_entry["conversion_factor"] is None and conv_key in definition:
                            try:
                                source_entry["conversion_factor"] = float(definition[conv_key])
                            except (TypeError, ValueError):
                                source_entry["conversion_factor"] = 1.0
        normalized = {}
        for clubb_name, entry in catalog.items():
            normalized[clubb_name] = {}
            for source_name in BENCHMARK_SOURCES:
                entries = tuple(
                    {"kind": item["kind"], "name": item["name"]}
                    for item in entry[source_name]["entries"]
                )
                if not entries:
                    continue
                conv = entry[source_name]["conversion_factor"]
                normalized[clubb_name][source_name] = {
                    "entries": entries,
                    "conversion_factor": 1.0 if conv is None else float(conv),
                }
        return normalized


def load_case_benchmark_catalog(case_data):
    benchmarks = (case_data or {}).get("benchmarks") or {}
    group_names = tuple(benchmarks.get("var_group_names") or [])
    case_name = (case_data or {}).get("name") or ""
    key = ("benchmark_catalog", case_name, group_names)

    def _build():
        return _build_benchmark_catalog(group_names)

    return _cached_lru(_CATALOG_CACHE, _CACHE_MAX_ENTRIES, key, _build)


def load_global_benchmark_catalog():
    return _cached_lru(
        _CATALOG_CACHE,
        _CACHE_MAX_ENTRIES,
        ("benchmark_catalog", "global"),
        lambda: _build_benchmark_catalog(_all_variable_group_names()),
    )


def _dataset_metadata(path):
    key = ("dataset_meta", _file_signature(path))

    def _build():
        ds = Dataset(path, "r")
        try:
            ds.set_auto_mask(False)
            var_info = {}
            vars_found = []
            dim_units = {}
            dim_names = tuple(ds.dimensions.keys())
            time_dim = _find_dim(dim_names, T_DIM_NAMES)
            time_values = None
            time_units = None
            if time_dim and time_dim in ds.variables:
                time_var = ds.variables[time_dim]
                time_values = np.asarray(time_var[:], dtype=float)
                time_units = getattr(time_var, "units", None)
                dim_units[time_dim] = getattr(time_var, "units", "") or ""
            for dim_name in dim_names:
                if dim_name in ds.variables and dim_name not in dim_units:
                    dim_units[dim_name] = getattr(ds.variables[dim_name], "units", "") or ""
            for name, var in ds.variables.items():
                if name in ds.dimensions:
                    continue
                dims = tuple(var.dimensions)
                vars_found.append(name)
                var_info[name] = {
                    "dims": dims,
                    "z_dim": _find_dim(dims, Z_DIM_NAMES),
                    "t_dim": _find_dim(dims, T_DIM_NAMES),
                    "units": getattr(var, "units", "") or "",
                    "long_name": getattr(var, "long_name", "") or "",
                }
            return {
                "vars": tuple(sorted(vars_found)),
                "var_info": var_info,
                "time_dim": time_dim,
                "time_values": _freeze(time_values) if time_values is not None else None,
                "time_units": time_units,
                "dim_units": dim_units,
            }
        finally:
            ds.close()

    return _cached_lru(_META_CACHE, _CACHE_MAX_ENTRIES, key, _build)


def _source_paths(case_data, source_name):
    return dict(((case_data or {}).get("benchmarks") or {}).get(source_name) or {})


def _ensure_profile_data(path, var_name, conversion_factor):
    key = ("profile_data", _file_signature(path), var_name, float(conversion_factor))

    def _build():
        meta = _dataset_metadata(path)
        info = meta["var_info"].get(var_name)
        if info is None or info.get("t_dim") is None or info.get("z_dim") is None:
            return None
        ds = Dataset(path, "r")
        try:
            ds.set_auto_mask(False)
            dims = info["dims"]
            t_dim = info["t_dim"]
            z_dim = info["z_dim"]
            var = ds.variables[var_name]
            slices = []
            kept_dims = []
            for dim in dims:
                if dim == t_dim or dim == z_dim:
                    slices.append(slice(None))
                    kept_dims.append(dim)
                    continue
                dim_size = len(ds.dimensions[dim])
                if dim_size > 1:
                    return None
                slices.append(0)
            data = np.asarray(var[tuple(slices)], dtype=float) * float(conversion_factor)
            if data.ndim != len(kept_dims):
                return None
            time_axis = kept_dims.index(t_dim)
            if time_axis != 0:
                data = np.moveaxis(data, time_axis, 0)
                kept_dims.insert(0, kept_dims.pop(time_axis))
            z_axis = kept_dims.index(z_dim)
            if z_axis != 1:
                data = np.moveaxis(data, z_axis, 1)
            z_var = ds.variables.get(z_dim)
            if z_var is not None:
                z_values = np.asarray(z_var[:], dtype=float)
                z_units = getattr(z_var, "units", "") or ""
            else:
                z_values = np.arange(data.shape[1], dtype=float)
                z_units = ""
            time_values = meta.get("time_values")
            if time_values is None:
                time_values = np.arange(data.shape[0], dtype=float)
                time_units = "seconds"
            else:
                time_units = meta.get("time_units")
            return {
                "profiles": np.asarray(data, dtype=float),
                "time_seconds": _time_values_seconds(time_values, time_units),
                "z_values": np.asarray(z_values, dtype=float),
                "units": info.get("units", "") or "",
                "long_name": info.get("long_name", "") or "",
                "z_units": z_units,
            }
        finally:
            ds.close()

    return _cached_lru(_DATA_CACHE, _CACHE_MAX_ENTRIES, key, lambda: _freeze(_build()))


def _ordered_path_items(source_paths, preferred_labels=None):
    preferred = [label for label in (preferred_labels or []) if label in source_paths]
    remaining = [label for label in source_paths if label not in preferred]
    return [(label, source_paths[label]) for label in preferred + remaining]


def _get_profile_data(source_paths, var_name, preferred_labels=None, conversion_factor=1.0):
    for label, path in _ordered_path_items(source_paths, preferred_labels=preferred_labels):
        meta = _dataset_metadata(path)
        info = meta["var_info"].get(var_name)
        if info is None or info.get("t_dim") is None or info.get("z_dim") is None:
            continue
        data = _ensure_profile_data(path, var_name, conversion_factor)
        if data:
            return data
    return None


def _first_profile(source_paths, candidates, preferred_labels=None, conversion_factor=1.0):
    for candidate in candidates:
        data = _get_profile_data(source_paths, candidate, preferred_labels=preferred_labels, conversion_factor=conversion_factor)
        if data:
            return data
    return None


def _compatible_profile_dict(data):
    return bool(data) and np.asarray(data["profiles"]).ndim == 2 and np.asarray(data["profiles"]).shape[0] > 0


def _same_grid(*datas):
    valid = [data for data in datas if _compatible_profile_dict(data)]
    if not valid:
        return False
    first = valid[0]
    first_profiles = np.asarray(first["profiles"], dtype=float)
    first_time = np.asarray(first["time_seconds"], dtype=float)
    first_z = np.asarray(first["z_values"], dtype=float)
    for data in valid[1:]:
        profiles = np.asarray(data["profiles"], dtype=float)
        if profiles.shape != first_profiles.shape:
            return False
        if not np.allclose(np.asarray(data["time_seconds"], dtype=float), first_time):
            return False
        if not np.allclose(np.asarray(data["z_values"], dtype=float), first_z):
            return False
    return True


def _required_same_grid(*datas):
    return all(_compatible_profile_dict(data) for data in datas) and _same_grid(*datas)


def _window_indices(data, clubb_window):
    benchmark_elapsed = _elapsed_seconds(data["time_seconds"])
    if clubb_window["mode"] == "point":
        index = int(np.argmin(np.abs(benchmark_elapsed - clubb_window["start"])))
        return np.array([index], dtype=int)
    low = min(clubb_window["start"], clubb_window["end"])
    high = max(clubb_window["start"], clubb_window["end"])
    matching = np.where((benchmark_elapsed >= low) & (benchmark_elapsed <= high))[0]
    if matching.size == 0:
        start_idx = int(np.argmin(np.abs(benchmark_elapsed - low)))
        end_idx = int(np.argmin(np.abs(benchmark_elapsed - high)))
        if end_idx < start_idx:
            start_idx, end_idx = end_idx, start_idx
        matching = np.arange(start_idx, end_idx + 1, dtype=int)
    return matching


def _reduce_profile_data(data, clubb_window):
    if not _compatible_profile_dict(data):
        return None
    indices = _window_indices(data, clubb_window)
    profiles = np.asarray(data["profiles"], dtype=float)
    if clubb_window["mode"] == "point":
        profile = profiles[int(indices[0]), :]
    else:
        profile = np.mean(profiles[indices, :], axis=0)
    return {
        "profile": np.asarray(profile, dtype=float),
        "z_values": np.asarray(data["z_values"], dtype=float),
        "units": data.get("units", "") or "",
        "long_name": data.get("long_name", "") or "",
        "z_units": data.get("z_units", "") or "",
    }


def _with_profile(base, profile, units=None, long_name=None):
    reduced = _reduce_profile_data(base, {"mode": "point", "start": 0.0, "end": 0.0})
    if reduced is None:
        return None
    reduced["profile"] = np.asarray(profile, dtype=float)
    if units is not None:
        reduced["units"] = units
    if long_name is not None:
        reduced["long_name"] = long_name
    return reduced


def _with_profiles(base, profiles, units=None, long_name=None):
    if not _compatible_profile_dict(base):
        return None
    result = dict(base)
    result["profiles"] = np.asarray(profiles, dtype=float)
    if units is not None:
        result["units"] = units
    if long_name is not None:
        result["long_name"] = long_name
    return _freeze(result)


def _optional_profiles(data):
    if not _compatible_profile_dict(data):
        return None
    profiles = np.asarray(data["profiles"], dtype=float)
    if np.all(np.isnan(profiles)):
        return None
    return profiles


def _pick_nonzero_output(output1, output2):
    if output1 is None:
        return output2
    if output2 is None:
        return output1
    out1 = np.asarray(output1, dtype=float)
    out2 = np.asarray(output2, dtype=float)
    out1_is_nan = np.isnan(np.nanmean(out1))
    out2_is_nan = np.isnan(np.nanmean(out2))
    out1_is_zero = np.all(np.isclose(out1, 0))
    out2_is_zero = np.all(np.isclose(out2, 0))
    if (out1_is_zero or out1_is_nan) and not out2_is_nan:
        return out2
    if (out2_is_zero or out2_is_nan) and not out1_is_nan:
        return out1
    return out1


def _apply_conversion(data, factor):
    factor = float(factor)
    if data is None or factor == 1.0:
        return data
    if "profile" in data:
        result = dict(data)
        result["profile"] = np.asarray(data["profile"], dtype=float) * factor
        return result
    if _compatible_profile_dict(data):
        return _with_profiles(data, np.asarray(data["profiles"], dtype=float) * factor)
    return data


def _zero_profile_like(data, units=None, long_name=None):
    reduced = _reduce_profile_data(data, {"mode": "point", "start": 0.0, "end": 0.0})
    if reduced is None:
        return None
    return {
        **reduced,
        "profile": np.zeros_like(np.asarray(reduced["profile"], dtype=float)),
        "units": reduced["units"] if units is None else units,
        "long_name": reduced["long_name"] if long_name is None else long_name,
    }


def _calc_sam_thlm(source_paths, clubb_window):
    thetal = _get_profile_data(source_paths, "THETAL")
    theta = _get_profile_data(source_paths, "THETA")
    tabs = _get_profile_data(source_paths, "TABS")
    qi = _get_profile_data(source_paths, "QI")
    if _required_same_grid(thetal, theta, tabs, qi):
        thetal_r = _reduce_profile_data(thetal, clubb_window)
        theta_r = _reduce_profile_data(theta, clubb_window)
        tabs_r = _reduce_profile_data(tabs, clubb_window)
        qi_r = _reduce_profile_data(qi, clubb_window)
        profile = np.asarray(thetal_r["profile"], dtype=float) + (2500.4 * (np.asarray(theta_r["profile"], dtype=float) / np.asarray(tabs_r["profile"], dtype=float)) * (np.asarray(qi_r["profile"], dtype=float) / 1000.0))
        return {**thetal_r, "profile": profile, "units": "K"}
    if _required_same_grid(thetal, theta):
        thetal_r = _reduce_profile_data(thetal, clubb_window)
        return {**thetal_r, "units": "K"}
    return None


def _calc_sam_rtm(source_paths, clubb_window):
    qt = _get_profile_data(source_paths, "QT")
    if not _compatible_profile_dict(qt):
        return None
    qi = _get_profile_data(source_paths, "QI")
    qt_r = _reduce_profile_data(qt, clubb_window)
    qt_profile = np.asarray(qt_r["profile"], dtype=float)
    if _required_same_grid(qt, qi):
        qi_r = _reduce_profile_data(qi, clubb_window)
        qi_profile = np.asarray(qi_r["profile"], dtype=float)
        profile = (qt_profile - qi_profile) / 1000.0
        if np.any(np.isnan(qi_profile)):
            profile = qt_profile / 1000.0
    else:
        profile = qt_profile / 1000.0
    return {**qt_r, "profile": profile, "units": "kg/kg"}


def _calc_skewness(source_paths, num_candidates, den_candidates, offset, clubb_window):
    numerator = _first_profile(source_paths, num_candidates)
    denominator = _first_profile(source_paths, den_candidates)
    if not _required_same_grid(numerator, denominator):
        return None
    numerator_r = _reduce_profile_data(numerator, clubb_window)
    denominator_r = _reduce_profile_data(denominator, clubb_window)
    profile = np.asarray(numerator_r["profile"], dtype=float) / (np.asarray(denominator_r["profile"], dtype=float) + offset) ** 1.5
    return {**numerator_r, "profile": profile}


def _calc_sam_wpthlp(source_paths, clubb_window):
    tlflux = _get_profile_data(source_paths, "TLFLUX")
    rho = _get_profile_data(source_paths, "RHO")
    if not _required_same_grid(tlflux, rho):
        return None
    tlflux_r = _reduce_profile_data(tlflux, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    profile = np.asarray(tlflux_r["profile"], dtype=float) / (np.asarray(rho_r["profile"], dtype=float) * 1004.0)
    sgs = _reduce_profile_data(_get_profile_data(source_paths, "WPTHLP_SGS"), clubb_window) if _compatible_profile_dict(_get_profile_data(source_paths, "WPTHLP_SGS")) else None
    if sgs is not None:
        profile = profile + np.asarray(sgs["profile"], dtype=float)
    return {**tlflux_r, "profile": profile}


def _calc_sam_wprtp(source_paths, clubb_window):
    qtflux = _get_profile_data(source_paths, "QTFLUX")
    rho = _get_profile_data(source_paths, "RHO")
    if not _required_same_grid(qtflux, rho):
        return None
    qtflux_r = _reduce_profile_data(qtflux, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    profile = np.asarray(qtflux_r["profile"], dtype=float) / (np.asarray(rho_r["profile"], dtype=float) * 2.5104e6)
    sgs_full = _get_profile_data(source_paths, "WPRTP_SGS")
    if _compatible_profile_dict(sgs_full):
        sgs = _reduce_profile_data(sgs_full, clubb_window)
        profile = profile + np.asarray(sgs["profile"], dtype=float)
    return {**qtflux_r, "profile": profile}


def _calc_sam_wpthvp(source_paths, clubb_window):
    tvflux = _get_profile_data(source_paths, "TVFLUX")
    rho = _get_profile_data(source_paths, "RHO")
    if not _required_same_grid(tvflux, rho):
        return None
    tvflux_r = _reduce_profile_data(tvflux, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    profile = np.asarray(tvflux_r["profile"], dtype=float) / (np.asarray(rho_r["profile"], dtype=float) * 1004.0)
    return {**tvflux_r, "profile": profile}


def _calc_sam_additive(source_paths, base_var, sgs_var, clubb_window):
    base = _get_profile_data(source_paths, base_var)
    if not _compatible_profile_dict(base):
        return None
    base_r = _reduce_profile_data(base, clubb_window)
    profile = np.asarray(base_r["profile"], dtype=float)
    sgs_full = _get_profile_data(source_paths, sgs_var)
    if _compatible_profile_dict(sgs_full):
        sgs_r = _reduce_profile_data(sgs_full, clubb_window)
        profile = profile + np.asarray(sgs_r["profile"], dtype=float)
    return {**base_r, "profile": profile}


def _calc_sam_rtp2(source_paths, clubb_window):
    qt2 = _get_profile_data(source_paths, "QT2")
    if not _compatible_profile_dict(qt2):
        return None
    qt2_r = _reduce_profile_data(qt2, clubb_window)
    profile = np.asarray(qt2_r["profile"], dtype=float) / 1.0e6
    sgs_full = _get_profile_data(source_paths, "RTP2_SGS")
    if _compatible_profile_dict(sgs_full):
        sgs_r = _reduce_profile_data(sgs_full, clubb_window)
        profile = profile + np.asarray(sgs_r["profile"], dtype=float)
    return {**qt2_r, "profile": profile, "units": "kg2/kg2"}


def _calc_sam_rtp3(source_paths, clubb_window):
    rc_coef = _get_profile_data(source_paths, "rc_coef_zm")
    rtprcp = _get_profile_data(source_paths, "rtprcp")
    if _required_same_grid(rc_coef, rtprcp):
        rc_r = _reduce_profile_data(rc_coef, clubb_window)
        rt_r = _reduce_profile_data(rtprcp, clubb_window)
        return {**rc_r, "profile": np.asarray(rc_r["profile"], dtype=float) * np.asarray(rt_r["profile"], dtype=float)}
    qcflux = _get_profile_data(source_paths, "QCFLUX")
    rho = _get_profile_data(source_paths, "RHO")
    pres = _get_profile_data(source_paths, "PRES")
    thetav = _get_profile_data(source_paths, "THETAV")
    if not _required_same_grid(qcflux, rho, pres, thetav):
        return None
    qcflux_r = _reduce_profile_data(qcflux, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    pres_r = _reduce_profile_data(pres, clubb_window)
    thetav_r = _reduce_profile_data(thetav, clubb_window)
    output = (np.asarray(qcflux_r["profile"], dtype=float) / (np.asarray(rho_r["profile"], dtype=float) * 2.5104e6)) * (
        2.5e6 / (1004.67 * ((np.asarray(pres_r["profile"], dtype=float) / 1000.0) ** (287.04 / 1004.67)))
        - 1.61 * np.asarray(thetav_r["profile"], dtype=float)
    )
    return {**qcflux_r, "profile": output}


def _calc_sam_tke(source_paths, clubb_window):
    u2 = _calc_sam_additive(source_paths, "U2", "UP2_SGS", clubb_window)
    v2 = _calc_sam_additive(source_paths, "V2", "VP2_SGS", clubb_window)
    w2 = _calc_sam_additive(source_paths, "W2", "WP2_SGS", clubb_window)
    if any(item is None for item in (u2, v2, w2)):
        return None
    profile = 0.5 * (np.asarray(u2["profile"], dtype=float) + np.asarray(v2["profile"], dtype=float) + np.asarray(w2["profile"], dtype=float))
    return {**u2, "profile": profile}


def _calc_sam_bv_freq_sqd(source_paths, clubb_window):
    thetav = _get_profile_data(source_paths, "THETAV")
    if not _compatible_profile_dict(thetav):
        return None
    thetav_r = _reduce_profile_data(thetav, clubb_window)
    z_values = np.asarray(thetav_r["z_values"], dtype=float)
    if z_values.size < 2:
        return None
    profile = np.asarray(thetav_r["profile"], dtype=float)
    tmp = np.diff(profile) / np.diff(z_values)
    deriv = np.zeros_like(profile)
    deriv[:-1] = tmp
    deriv[-1] = tmp[-1]
    output = 9.81 / profile * deriv
    return {**thetav_r, "profile": output, "units": "1/s^2"}


def _calc_sam_ncm(source_paths, clubb_window):
    rho = _get_profile_data(source_paths, "RHO")
    cld = _get_profile_data(source_paths, "CLD")
    if not _required_same_grid(rho, cld):
        return None
    rho_r = _reduce_profile_data(rho, clubb_window)
    cld_r = _reduce_profile_data(cld, clubb_window)
    nc = _get_profile_data(source_paths, "NC")
    gcssnc = _get_profile_data(source_paths, "GCSSNC")
    output1 = None
    output2 = None
    if _required_same_grid(nc, rho, cld):
        nc_r = _reduce_profile_data(nc, clubb_window)
        output1 = np.asarray(cld_r["profile"], dtype=float) * (np.asarray(nc_r["profile"], dtype=float) * 1.0e6 / np.asarray(rho_r["profile"], dtype=float))
    if _required_same_grid(gcssnc, rho):
        gcssnc_r = _reduce_profile_data(gcssnc, clubb_window)
        output2 = np.asarray(gcssnc_r["profile"], dtype=float) * 1.0e6 / np.asarray(rho_r["profile"], dtype=float)
    chosen = _pick_nonzero_output(output1, output2)
    if chosen is None:
        return _zero_profile_like(rho, units="#/kg", long_name="Ncm")
    return {**rho_r, "profile": chosen}


def _calc_sam_nc_in_cloud(source_paths, clubb_window):
    nc = _get_profile_data(source_paths, "NC")
    rho = _get_profile_data(source_paths, "RHO")
    if not _compatible_profile_dict(rho):
        return None
    if not _required_same_grid(nc, rho):
        return _zero_profile_like(rho, units="#/kg", long_name="Nc_in_cloud")
    nc_r = _reduce_profile_data(nc, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    output = np.asarray(nc_r["profile"], dtype=float) * 1.0e6 / np.asarray(rho_r["profile"], dtype=float)
    return {**nc_r, "profile": output}


def _calc_sam_nrm(source_paths, clubb_window):
    nr = _first_profile(source_paths, ["NR", "CONP"])
    rho = _get_profile_data(source_paths, "RHO")
    if not _compatible_profile_dict(rho):
        return None
    if not _required_same_grid(nr, rho):
        return _zero_profile_like(rho, units="#/kg", long_name="Nrm")
    nr_r = _reduce_profile_data(nr, clubb_window)
    rho_r = _reduce_profile_data(rho, clubb_window)
    output = np.asarray(nr_r["profile"], dtype=float) * 1.0e6 / np.asarray(rho_r["profile"], dtype=float)
    return {**nr_r, "profile": output}


def _calc_coamps_upwp(source_paths, clubb_window):
    wpup = _get_profile_data(source_paths, "wpup", preferred_labels=["sw"])
    wpup_sgs = _get_profile_data(source_paths, "wpup_sgs", preferred_labels=["sw"])
    if not _required_same_grid(wpup, wpup_sgs):
        return None
    wpup_r = _reduce_profile_data(wpup, clubb_window)
    wpup_sgs_r = _reduce_profile_data(wpup_sgs, clubb_window)
    return {**wpup_r, "profile": np.asarray(wpup_r["profile"], dtype=float) + np.asarray(wpup_sgs_r["profile"], dtype=float)}


def _calc_coamps_vpwp(source_paths, clubb_window):
    wpvp = _get_profile_data(source_paths, "wpvp", preferred_labels=["sw"])
    wpvp_sgs = _get_profile_data(source_paths, "wpvp_sgs", preferred_labels=["sw"])
    if not _required_same_grid(wpvp, wpvp_sgs):
        return None
    wpvp_r = _reduce_profile_data(wpvp, clubb_window)
    wpvp_sgs_r = _reduce_profile_data(wpvp_sgs, clubb_window)
    return {**wpvp_r, "profile": np.asarray(wpvp_r["profile"], dtype=float) + np.asarray(wpvp_sgs_r["profile"], dtype=float)}


_DERIVED_PROFILE_REGISTRY = {
    ("sam", "getThlmSamCalc"): _calc_sam_thlm,
    ("sam", "getRtmSamCalc"): _calc_sam_rtm,
    ("sam", "getSkwZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["WP3", "W3", "wp3"], ["WP2", "W2", "wp2"], 1.6e-3, window),
    ("coamps", "getSkwZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["WP3", "W3", "wp3"], ["WP2", "W2", "wp2"], 1.6e-3, window),
    ("sam", "getSkrtZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["RTP3", "qtp3", "rtp3"], ["RTP2", "qtp2", "rtp2", "rlp2"], 4e-16, window),
    ("coamps", "getSkrtZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["RTP3", "qtp3", "rtp3"], ["RTP2", "qtp2", "rtp2", "rlp2"], 4e-16, window),
    ("sam", "getSkthlZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["THLP3", "thlp3"], ["THLP2", "thlp2"], 4e-4, window),
    ("coamps", "getSkthlZtLesCalc"): lambda paths, window: _calc_skewness(paths, ["THLP3", "thlp3"], ["THLP2", "thlp2"], 4e-4, window),
    ("sam", "getWpthlpSamCalc"): _calc_sam_wpthlp,
    ("sam", "getWprtpSamCalc"): _calc_sam_wprtp,
    ("sam", "getWpthvpSamCalc"): _calc_sam_wpthvp,
    ("sam", "get_wp2_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "W2", "WP2_SGS", window),
    ("sam", "get_wp3_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "W3", "WP3_SGS", window),
    ("sam", "get_thlp2_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "TL2", "THLP2_SGS", window),
    ("sam", "getRtp2SamCalc"): _calc_sam_rtp2,
    ("sam", "getRtp3SamCalc"): _calc_sam_rtp3,
    ("sam", "get_upwp_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "UW", "UPWP_SGS", window),
    ("sam", "get_vpwp_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "VW", "VPWP_SGS", window),
    ("sam", "get_up2_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "U2", "UP2_SGS", window),
    ("sam", "get_vp2_sam_calc"): lambda paths, window: _calc_sam_additive(paths, "V2", "VP2_SGS", window),
    ("sam", "get_tke_sam_calc"): _calc_sam_tke,
    ("sam", "getBVSqdSamCalc"): _calc_sam_bv_freq_sqd,
    ("sam", "getNcmSamLine"): _calc_sam_ncm,
    ("sam", "getNcInCloudSamLine"): _calc_sam_nc_in_cloud,
    ("sam", "getNrmSamLine"): _calc_sam_nrm,
    ("coamps", "getUwCoampsData"): _calc_coamps_upwp,
    ("coamps", "getVwCoampsData"): _calc_coamps_vpwp,
}


def _extract_from_catalog_entry(case_data, source_name, catalog_entry, clubb_window):
    source_paths = _source_paths(case_data, source_name)
    if not source_paths:
        return None
    if catalog_entry["kind"] == "calc":
        builder = _DERIVED_PROFILE_REGISTRY.get((source_name, catalog_entry["name"]))
        if builder is None:
            return None
        return builder(source_paths, clubb_window)
    if catalog_entry["kind"] == "var":
        return _reduce_profile_data(_get_profile_data(source_paths, catalog_entry["name"]), clubb_window)
    return None


def _resolve_clubb_elapsed(case_data, time_mode, time_range, time_point):
    time_seconds = np.asarray((case_data or {}).get("time_seconds") or [], dtype=float)
    time_len = max(int((case_data or {}).get("time_len") or 1), 1)
    if time_seconds.size > 0:
        elapsed = _elapsed_seconds(time_seconds)
        if time_mode == "point":
            idx = _slider_value_to_index(time_point, len(elapsed))
            return {"mode": "point", "start": float(elapsed[idx]), "end": float(elapsed[idx])}
        start_idx, end_idx = _slider_range_to_indices(time_range, len(elapsed))
        return {"mode": "range", "start": float(elapsed[start_idx]), "end": float(elapsed[end_idx])}
    if time_mode == "point":
        idx = _slider_value_to_index(time_point, time_len)
        return {"mode": "point", "start": float(idx), "end": float(idx)}
    start_idx, end_idx = _slider_range_to_indices(time_range, time_len)
    return {"mode": "range", "start": float(start_idx), "end": float(end_idx)}


def extract_benchmark_profile(case_data, source_name, clubb_var, time_mode, time_range, time_point):
    catalog = load_case_benchmark_catalog(case_data)
    source_catalog = (catalog.get(clubb_var) or {}).get(source_name)
    if not source_catalog:
        source_catalog = (load_global_benchmark_catalog().get(clubb_var) or {}).get(source_name)
    if not source_catalog:
        return None
    clubb_window = _resolve_clubb_elapsed(case_data, time_mode, time_range, time_point)
    data = None
    for catalog_entry in source_catalog["entries"]:
        candidate = _extract_from_catalog_entry(case_data, source_name, catalog_entry, clubb_window)
        candidate = _apply_conversion(candidate, source_catalog["conversion_factor"])
        if candidate is not None and "profile" in candidate:
            data = candidate
            break
    if data is None:
        return None
    return {
        "z_values": np.asarray(data["z_values"], dtype=float),
        "profile": np.asarray(data["profile"], dtype=float),
        "label": _profile_label(source_name),
        "line": _benchmark_line_style(source_name),
        "units": data["units"],
        "long_name": data["long_name"],
        "z_units": data["z_units"],
    }
