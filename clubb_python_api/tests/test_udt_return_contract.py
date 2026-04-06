"""Guardrails for Python UDT return semantics in clubb_python wrappers."""

import ast
import importlib
import inspect
import sys
import textwrap
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.derived_types.config_flags import ConfigFlags
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.pdf_params import implicit_coefs_terms, pdf_parameter
from clubb_python.derived_types.sclr_idx import SclrIdx


UDT_TYPES = (
    Grid,
    ErrInfo,
    NuVertResDep,
    ConfigFlags,
    SclrIdx,
    pdf_parameter,
    implicit_coefs_terms,
)

DIRECT_GETTER_TO_KIND = {
    "get_fortran_grid": "grid",
    "get_fortran_err_info": "err_info",
    "get_fortran_nu_vert_res_dep": "nu_vert_res_dep",
    "get_fortran_sclr_idx": "sclr_idx",
    "get_fortran_pdf_params": "pdf_params",
    "get_fortran_pdf_params_zm": "pdf_params_zm",
    "get_fortran_implicit_coefs": "implicit_coefs_terms",
}

UDT_ARG_RETURN_MANIFEST = {
    "advance_clubb_core.advance_clubb_core": {("pdf_params", "pdf_params_zm", "implicit_coefs_terms", "err_info")},
    "advance_windm_edsclrm.advance_windm_edsclrm": {("err_info",)},
    "advance_wp2_wp3.advance_wp2_wp3": {("err_info",)},
    "advance_xm_wpxp.advance_xm_wpxp": {("err_info",)},
    "advance_xp2_xpyp.advance_xp2_xpyp": {("err_info",)},
    "corr_varnce_module.assert_corr_symmetric": {("err_info",)},
    "grid_class.read_grid_heights": {("err_info",)},
    "grid_class.setup_grid_heights": {("grid", "err_info")},
    "lapack_wrap.lapack_band_solve": {("err_info",)},
    "lapack_wrap.lapack_band_solvex": {("err_info",)},
    "lapack_wrap.lapack_tridiag_solve": {("err_info",)},
    "lapack_wrap.lapack_tridiag_solvex": {("err_info",)},
    "matrix_solver_wrapper.band_solve": {("err_info",)},
    "matrix_solver_wrapper.tridiag_solve": {("err_info",)},
    "mixing_length.calc_lscale_directly": {("err_info",)},
    "mixing_length.diagnose_lscale_from_tau": {("err_info",)},
    "numerical_check.check_clubb_settings": {("err_info",)},
    "numerical_check.length_check": {("err_info",)},
    "numerical_check.parameterization_check": {("err_info",)},
    "numerical_check.pdf_closure_check": {("err_info",)},
    "numerical_check.sfc_varnce_check": {("err_info",)},
    "parameters_tunable.calc_derrived_params": {("nu_vert_res_dep",)},
    "parameters_tunable.check_parameters": {("err_info",)},
    "pdf_closure.pdf_closure_driver": {("implicit_coefs_terms", "pdf_params", "pdf_params_zm", "err_info")},
    "pdf_parameter_module.zero_pdf_params": {
        ("pdf_params",),
        ("pdf_params_zm",),
        ("pdf_params", "pdf_params_zm"),
    },
    "precipitation_fraction.precip_fraction": {("err_info",)},
    "grid_class.setup_grid": {("grid", "err_info")},
    "sfc_varnce_module.calc_sfc_varnce": {("err_info",)},
    "stats_netcdf.finalize_stats": {("err_info",)},
    "stats_netcdf.init_stats": {("err_info",)},
    "stats_netcdf.stats_end_timestep": {("err_info",)},
}

NOARG_UDT_RETURN_MANIFEST = {
    "pdf_parameter_module.init_pdf_implicit": {("implicit_coefs_terms",)},
    "pdf_parameter_module.init_pdf_params": {("pdf_params",)},
    "pdf_parameter_module.init_pdf_params_zm": {("pdf_params_zm",)},
}

FULL_RETURN_PATTERN_MANIFEST = {
    "advance_clubb_core.advance_clubb_core": {
        ("*", "pdf_params", "pdf_params_zm", "implicit_coefs_terms", "err_info", "*"),
    },
    "advance_windm_edsclrm.advance_windm_edsclrm": {
        ("*", "err_info"),
    },
    "advance_wp2_wp3.advance_wp2_wp3": {
        ("*", "err_info"),
    },
    "advance_xm_wpxp.advance_xm_wpxp": {
        ("*", "err_info"),
    },
    "advance_xp2_xpyp.advance_xp2_xpyp": {
        ("*", "err_info"),
    },
    "corr_varnce_module.assert_corr_symmetric": {
        ("err_info",),
    },
    "grid_class.read_grid_heights": {
        ("err_info", "value", "value"),
    },
    "grid_class.setup_grid_heights": {
        ("grid", "err_info"),
    },
    "lapack_wrap.lapack_band_solve": {
        ("err_info", "value"),
    },
    "lapack_wrap.lapack_band_solvex": {
        ("err_info", "*"),
    },
    "lapack_wrap.lapack_tridiag_solve": {
        ("err_info", "value"),
    },
    "lapack_wrap.lapack_tridiag_solvex": {
        ("err_info", "*"),
    },
    "matrix_solver_wrapper.band_solve": {
        ("err_info", "value", "value"),
    },
    "matrix_solver_wrapper.tridiag_solve": {
        ("err_info", "value", "value"),
    },
    "mixing_length.calc_lscale_directly": {
        ("err_info", "*"),
    },
    "mixing_length.diagnose_lscale_from_tau": {
        ("err_info", "*"),
    },
    "numerical_check.check_clubb_settings": {
        ("err_info",),
    },
    "numerical_check.length_check": {
        ("err_info",),
    },
    "numerical_check.parameterization_check": {
        ("err_info",),
    },
    "numerical_check.pdf_closure_check": {
        ("err_info",),
    },
    "numerical_check.sfc_varnce_check": {
        ("err_info",),
    },
    "parameters_tunable.calc_derrived_params": {
        ("nu_vert_res_dep", "value", "value"),
    },
    "parameters_tunable.check_parameters": {
        ("err_info",),
    },
    "pdf_closure.pdf_closure_driver": {
        ("value", "implicit_coefs_terms", "pdf_params", "pdf_params_zm", "err_info", "*"),
    },
    "pdf_parameter_module.init_pdf_implicit": {
        ("implicit_coefs_terms",),
    },
    "pdf_parameter_module.init_pdf_params": {
        ("pdf_params",),
    },
    "pdf_parameter_module.init_pdf_params_zm": {
        ("pdf_params_zm",),
    },
    "prescribe_forcings.prescribe_forcings": {
        ("*", "err_info"),
    },
    "precipitation_fraction.precip_fraction": {
        ("err_info", "value", "value", "value", "value"),
    },
    "grid_class.setup_grid": {
        ("grid", "err_info"),
    },
    "sfc_varnce_module.calc_sfc_varnce": {
        ("*", "err_info"),
    },
    "stats_netcdf.finalize_stats": {
        ("err_info",),
    },
    "stats_netcdf.init_stats": {
        ("err_info",),
    },
    "stats_netcdf.stats_end_timestep": {
        ("err_info",),
    },
}

SOURCE_ORDER_MANIFEST = {
    "advance_clubb_core.advance_clubb_core": {
        "source": "src/CLUBB_core/advance_clubb_core_module.F90::advance_clubb_core",
        "hidden_removed": (
            "gr",
            "stored_sclr_idx",
            "stored_nu_vert_res_dep",
            "stored_config_flags",
            "stored_stats",
            "stored_pdf_params",
            "stored_pdf_params_zm",
            "stored_pdf_implicit_coefs_terms",
            "stored_err_info",
        ),
        "return_prefix_count": 51,
        "udt_return_order": ("pdf_params", "pdf_params_zm", "implicit_coefs_terms", "err_info"),
    },
    "pdf_closure.pdf_closure_driver": {
        "source": "src/CLUBB_core/pdf_closure_module.F90::pdf_closure_driver",
        "hidden_removed": ("gr", "stats", "pdf_implicit_coefs_terms", "pdf_params", "pdf_params_zm", "err_info"),
        "udt_return_order": ("implicit_coefs_terms", "pdf_params", "pdf_params_zm", "err_info"),
    },
}


def _call_name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return None


def _annotation_contains_udt(annotation) -> bool:
    if annotation is inspect._empty:
        return False
    if isinstance(annotation, str):
        return annotation in {cls.__name__ for cls in UDT_TYPES}
    if annotation in UDT_TYPES:
        return True
    args = getattr(annotation, "__args__", None)
    if args:
        return any(_annotation_contains_udt(arg) for arg in args)
    return False


def _extract_return_udt_sequences(func) -> set[tuple[str, ...]]:
    source = textwrap.dedent(inspect.getsource(func))
    tree = ast.parse(source)
    sequences: set[tuple[str, ...]] = set()

    def extract_kinds(expr) -> list[str]:
        if isinstance(expr, ast.Call):
            call_name = _call_name(expr.func)
            if call_name in DIRECT_GETTER_TO_KIND:
                return [DIRECT_GETTER_TO_KIND[call_name]]
            return []
        if isinstance(expr, ast.Tuple):
            kinds = []
            for elt in expr.elts:
                if isinstance(elt, ast.Starred):
                    kinds.extend(extract_kinds(elt.value))
                else:
                    kinds.extend(extract_kinds(elt))
            return kinds
        return []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Return) or node.value is None:
            continue
        udt_kinds = tuple(extract_kinds(node.value))
        if udt_kinds:
            sequences.add(udt_kinds)

    return sequences


def _extract_return_patterns(func) -> set[tuple[str, ...]]:
    source = textwrap.dedent(inspect.getsource(func))
    tree = ast.parse(source)
    patterns: set[tuple[str, ...]] = set()

    def extract_pattern(expr) -> list[str]:
        if isinstance(expr, ast.Call):
            call_name = _call_name(expr.func)
            if call_name in DIRECT_GETTER_TO_KIND:
                return [DIRECT_GETTER_TO_KIND[call_name]]
            return ["value"]
        if isinstance(expr, ast.Tuple):
            pattern = []
            for elt in expr.elts:
                if isinstance(elt, ast.Starred):
                    pattern.append("*")
                else:
                    pattern.extend(extract_pattern(elt))
            return pattern
        return ["value"]

    for node in ast.walk(tree):
        if not isinstance(node, ast.Return) or node.value is None:
            continue
        patterns.add(tuple(extract_pattern(node.value)))

    return patterns


def _extract_return_slice_markers(func) -> set[tuple[tuple[str, int], ...]]:
    source = textwrap.dedent(inspect.getsource(func))
    tree = ast.parse(source)
    markers: set[tuple[tuple[str, int], ...]] = set()

    def extract(expr):
        if not isinstance(expr, ast.Tuple):
            return tuple()
        out = []
        for elt in expr.elts:
            if isinstance(elt, ast.Starred) and isinstance(elt.value, ast.Subscript):
                sub = elt.value
                if isinstance(sub.value, ast.Name) and sub.value.id == "result" and isinstance(sub.slice, ast.Slice):
                    lower = sub.slice.lower.value if isinstance(sub.slice.lower, ast.Constant) else None
                    upper = sub.slice.upper.value if isinstance(sub.slice.upper, ast.Constant) else None
                    if lower is None and upper is not None:
                        out.append(("prefix", int(upper)))
                    elif lower is not None and upper is None:
                        out.append(("suffix", int(lower)))
            elif isinstance(elt, ast.Subscript):
                if isinstance(elt.value, ast.Name) and elt.value.id == "result" and isinstance(elt.slice, ast.Constant):
                    out.append(("index", int(elt.slice.value)))
        return tuple(out)

    for node in ast.walk(tree):
        if isinstance(node, ast.Return) and node.value is not None:
            markers.add(extract(node.value))

    return markers


def _load_clubb_core_modules():
    pkg_dir = Path(__file__).resolve().parent.parent / "clubb_python" / "CLUBB_core"
    modules = {}
    for path in sorted(pkg_dir.glob("*.py")):
        if path.name == "__init__.py":
            continue
        stem = path.stem
        modules[stem] = importlib.import_module(f"clubb_python.CLUBB_core.{stem}")
    return modules


def _discover_wrappers_with_udt_args():
    discovered = {}
    for stem, module in _load_clubb_core_modules().items():
        for name, func in inspect.getmembers(module, inspect.isfunction):
            if func.__module__ != module.__name__:
                continue
            if any(_annotation_contains_udt(param.annotation) for param in inspect.signature(func).parameters.values()):
                discovered[f"{stem}.{name}"] = func
    return discovered


def _load_function(qualified_name: str):
    module_name, func_name = qualified_name.split(".", 1)
    try:
        module = importlib.import_module(f"clubb_python.CLUBB_core.{module_name}")
    except ModuleNotFoundError:
        module = importlib.import_module(f"clubb_python.{module_name}")
    return getattr(module, func_name)


def test_udt_arg_wrapper_return_manifest_matches_source():
    """Wrappers with UDT args must append refreshed UDTs in Fortran inout/out order."""
    discovered = _discover_wrappers_with_udt_args()
    assert discovered, "Expected to discover at least one wrapper with Python UDT arguments."

    unexpected_manifest_keys = set(UDT_ARG_RETURN_MANIFEST) - set(discovered)
    assert not unexpected_manifest_keys, (
        "Manifest includes wrappers not discovered as having UDT args: "
        + ", ".join(sorted(unexpected_manifest_keys))
    )

    for qualified_name, func in sorted(discovered.items()):
        actual = _extract_return_udt_sequences(func)
        expected = UDT_ARG_RETURN_MANIFEST.get(qualified_name, set())
        assert actual == expected, (
            f"{qualified_name} has UDT return sequences {sorted(actual)} "
            f"but expected {sorted(expected)}"
        )


def test_noarg_udt_return_wrappers_match_source():
    """Selected no-UDT-arg wrappers that create UDT state must return it explicitly."""
    for qualified_name, expected in sorted(NOARG_UDT_RETURN_MANIFEST.items()):
        actual = _extract_return_udt_sequences(_load_function(qualified_name))
        assert actual == expected, (
            f"{qualified_name} has UDT return sequences {sorted(actual)} "
            f"but expected {sorted(expected)}"
        )


def test_full_return_patterns_match_source():
    """Mixed wrappers must preserve full Fortran inout/out ordering, not just UDT sub-order."""
    for qualified_name, expected in sorted(FULL_RETURN_PATTERN_MANIFEST.items()):
        actual = _extract_return_patterns(_load_function(qualified_name))
        assert actual == expected, (
            f"{qualified_name} has return patterns {sorted(actual)} "
            f"but expected {sorted(expected)}"
        )


def test_source_order_manifest_matches_wrapper_slices():
    """Selected wrappers must splice numeric results around hidden UDTs at source-order positions."""
    for qualified_name, expected in sorted(SOURCE_ORDER_MANIFEST.items()):
        actual_udts = _extract_return_udt_sequences(_load_function(qualified_name))
        assert actual_udts == {expected["udt_return_order"]}

        actual_markers = _extract_return_slice_markers(_load_function(qualified_name))
        if "return_prefix_count" in expected:
            assert actual_markers == {
                (("prefix", expected["return_prefix_count"]), ("suffix", expected["return_prefix_count"]))
            }, (
                f"{qualified_name} has slice markers {sorted(actual_markers)} "
                f"but expected prefix/suffix split at {expected['return_prefix_count']}"
            )
        else:
            assert actual_markers == {(("index", 0), ("suffix", 1))}, (
                f"{qualified_name} has slice markers {sorted(actual_markers)} "
                "but expected an explicit leading result[0] followed by *result[1:]"
            )
