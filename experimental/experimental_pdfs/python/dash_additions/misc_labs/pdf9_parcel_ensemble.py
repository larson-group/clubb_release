"""Interactive diagnosis of PDF-9 parcel crossings as a transport Gaussian."""

from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any
from types import SimpleNamespace

from dash import Input, Output, State, ctx, dcc, html, no_update
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

if __package__ and __package__.startswith("dash_app."):
    from dash_app.shared.bivariate_heatmap import (
        BIVARIATE_LEVELS,
        signed_bivariate_codes,
        signed_bivariate_reference_colorscale,
    )
    from dash_app.shared.reporting import empty_state, fact_grid, report_header, report_section
else:
    from shared.bivariate_heatmap import (
        BIVARIATE_LEVELS,
        signed_bivariate_codes,
        signed_bivariate_reference_colorscale,
    )
    from shared.reporting import empty_state, fact_grid, report_header, report_section

from utilities.generate_pdf9_parcel_ensemble_diagnostics import (
    REFERENCE_NAMES,
    SCHEMA_VERSION,
)
from utilities.les_chi_moments import derive_chi_moments
from utilities.sam_3d_reference import (
    diagnose_cloud,
    load_snapshot,
    load_snapshot_pressure_profile,
)

from ..generation import job_status, start_pdf9_generation
from ..registry import ReportSpec, register_report


ROOT = Path(__file__).resolve().parents[3]
CACHE_DIR = ROOT / "output_pdf_bakeoff/pdf9_parcel_ensemble_lab"
MANIFEST_PATH = CACHE_DIR / "manifest.json"

CASE_ID = "notes-pdf9-case"
TIME_ID = "notes-pdf9-time"
HEIGHT_ID = "notes-pdf9-height"
OVERLAY_VARIABLE_ID = "notes-pdf9-overlay-variable"
TIME_LABEL_ID = "notes-pdf9-time-label"
HEIGHT_LABEL_ID = "notes-pdf9-height-label"
FACTS_ID = "notes-pdf9-facts"
OVERLAY_ID = "notes-pdf9-overlay"
COMPONENT_WEIGHTS_ID = "notes-pdf9-component-weights"
TRACKS_ID = "notes-pdf9-tracks"
REACH_ID = "notes-pdf9-reach"
SOURCES_ID = "notes-pdf9-sources"
ERROR_ID = "notes-pdf9-errors"
ANOMALY_ID = "notes-pdf9-anomaly"
BAKEOFF_ID = "notes-pdf9-bakeoff"
OCCUPANCY_COMPARISON_ID = "notes-pdf9-occupancy-comparison"
CLOSURE_ID = "notes-pdf9-closure"
CONTINUITY_ID = "notes-pdf9-continuity"
RECURSIVE_ID = "notes-pdf9-recursive"
ENTRAINMENT_ID = "notes-pdf9-entrainment"
CAP_ID = "notes-pdf9-downward-arrival-cap"
REFRESH_ID = "notes-pdf9-refresh-cache"
CACHE_STATUS_ID = "notes-pdf9-cache-status"
CACHE_REVISION_ID = "notes-pdf9-cache-revision"
GENERATION_INTERVAL_ID = "notes-pdf9-generation-interval"

WEIGHTINGS = {
    "uniform": "Occupancy (equal per crossing)",
    "speed": "Crossing speed w",
    "kinetic": "Kinetic proxy w²",
    "moist_flux": "Positive moist-anomaly flux",
    "buoyancy_flux": "Positive buoyancy flux",
}
OFFLINE_REFERENCE_WEIGHTING = "speed"
CLOSURE_PACKET_FRACTION = 0.08
OVERLAY_VARIABLES = {
    "rt": {"label": "w–rₜ", "axis": "Total water rₜ [g/kg]"},
    "chi": {"label": "w–χ", "axis": "Extended cloud water χ [g/kg]"},
    "rc": {"label": "w–r_c", "axis": "Cloud water r_c [g/kg]"},
}
COLORS = {
    "candidate": "#34d399",
    "sam": "#fbbf24",
    "mean": "#94a3b8",
    "up": "#22d3ee",
    "down": "#fb7185",
    "purple": "#a855f7",
    "recursive": "#bef264",
    "g1": "#38bdf8",
    "g2": "#f472b6",
}


def _empty(title: str, message: str = "No finite data at this selection"):
    figure = go.Figure()
    figure.add_annotation(text=message, showarrow=False, font={"color": "#94a3b8"})
    figure.update_layout(template="plotly_dark", title=title, height=420)
    return figure


def _finite_float(value: Any, default: float = 0.0) -> float:
    """Return one cached scalar without leaking masked/NaN boundary values."""
    if np.ma.is_masked(value):
        return default
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return default
    return scalar if np.isfinite(scalar) else default


def _finite_int(value: Any, default: int = 0) -> int:
    """Return an integer-valued diagnostic count, defaulting when unavailable."""
    scalar = _finite_float(value, math.nan)
    return int(scalar) if np.isfinite(scalar) else default


def _distance_label(value: float) -> str:
    return f"{value:.0f} m" if np.isfinite(value) else "unavailable"


@lru_cache(maxsize=1)
def _manifest() -> dict[str, Any]:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"No PDF-9 parcel cache exists at {MANIFEST_PATH}. Run: "
            "python utilities/generate_pdf9_parcel_ensemble_diagnostics.py"
        )
    payload = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if int(payload.get("schema_version", -1)) != SCHEMA_VERSION:
        raise ValueError("The PDF-9 parcel cache is stale; regenerate it")
    if not payload.get("cases"):
        raise ValueError("The PDF-9 parcel cache contains no cases")
    return payload


def _validate_cache_files(manifest: dict[str, Any]) -> None:
    """Fail early with an actionable message when a copied cache is incomplete."""

    missing = [
        str(CACHE_DIR / str(case.get("data_file", "<missing data_file>")))
        for case in manifest["cases"]
        if not case.get("data_file")
        or not (CACHE_DIR / str(case["data_file"])).is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "PDF-9 cache is incomplete; missing "
            + ", ".join(missing)
            + ". Click Generate / refresh PDF-9 data."
        )


@lru_cache(maxsize=4)
def _load_cube(path: str, mtime_ns: int) -> dict[str, np.ndarray]:
    del mtime_ns
    with np.load(path, allow_pickle=False) as source:
        return {name: np.asarray(source[name]).copy() for name in source.files}


def _case(name: str):
    manifest = _manifest()
    meta = next((item for item in manifest["cases"] if item["name"] == name), None)
    if meta is None:
        meta = manifest["cases"][0]
    path = CACHE_DIR / meta["data_file"]
    return meta, _load_cube(str(path), path.stat().st_mtime_ns)


def _clear_data_caches() -> None:
    """Forget every object derived from the replaceable on-disk cache."""
    _manifest.cache_clear()
    _load_cube.cache_clear()
    _raw_snapshot.cache_clear()
    _continuity_fields.cache_clear()
    _whole_case_center_metrics.cache_clear()
    load_snapshot_pressure_profile.cache_clear()


@lru_cache(maxsize=32)
def _raw_snapshot(run_dir: str, seconds: int, height_m: float):
    return load_snapshot(Path(run_dir), int(seconds), float(height_m))


def _weights(cube: dict[str, np.ndarray], ti: int, di: int, scheme: str):
    valid = cube["parcel_status"][ti, di] > 0
    w = np.maximum(cube["parcel_w"][ti, di].astype(float), 0.0)
    rt = cube["parcel_rt"][ti, di].astype(float)
    buoyancy = cube["parcel_buoyancy"][ti, di].astype(float)
    if scheme == "uniform":
        # Occupancy hypothesis: every successful parcel crossing represents
        # one equally occupied transport path.  The same weights are consumed
        # by _ensemble for all three means and the full covariance matrix.
        weight = np.ones_like(w)
    elif scheme == "kinetic":
        weight = w**2
    elif scheme == "moist_flux":
        weight = w * np.maximum(rt - float(cube["mean_rt"][ti, di]), 0.0)
    elif scheme == "buoyancy_flux":
        weight = w * np.maximum(buoyancy, 0.0)
    else:
        weight = w
    return np.where(valid & np.isfinite(weight), np.maximum(weight, 0.0), 0.0)


def _ensemble(cube, ti: int, di: int, scheme: str):
    weight = _weights(cube, ti, di, scheme)
    total = float(np.sum(weight))
    values = np.column_stack(
        (
            cube["parcel_w"][ti, di],
            cube["parcel_rt"][ti, di],
            cube["parcel_thl"][ti, di],
        )
    ).astype(float)
    finite = (weight > 0.0) & np.all(np.isfinite(values), axis=1)
    total = float(np.sum(weight[finite]))
    if total <= 1.0e-24:
        return {
            "available": False,
            "weight": weight,
            "mean": np.full(3, np.nan),
            "covariance": np.full((3, 3), np.nan),
            "neff": 0.0,
            "count": 0,
        }
    normalized = weight[finite] / total
    selected = values[finite]
    mean = np.sum(normalized[:, None] * selected, axis=0)
    anomaly = selected - mean
    covariance = np.einsum("n,ni,nj->ij", normalized, anomaly, anomaly)
    return {
        "available": True,
        "weight": weight,
        "mean": mean,
        "covariance": covariance,
        "neff": total**2 / max(float(np.sum(weight[finite] ** 2)), 1.0e-30),
        "count": int(np.count_nonzero(finite)),
    }


def _downward_ensemble(cube, ti: int, di: int):
    """Reduce the signed-w downward ledger with sinking-speed weights."""
    valid = cube["down_parcel_status"][ti, di] > 0
    values = np.column_stack(
        (
            cube["down_parcel_w"][ti, di],
            cube["down_parcel_rt"][ti, di],
            cube["down_parcel_thl"][ti, di],
        )
    ).astype(float)
    weight = np.where(valid, np.maximum(-values[:, 0], 0.0), 0.0)
    finite = (weight > 0.0) & np.all(np.isfinite(values), axis=1)
    total = float(np.sum(weight[finite]))
    if total <= 1.0e-24:
        return {
            "available": False,
            "weight": weight,
            "mean": np.full(3, np.nan),
            "covariance": np.full((3, 3), np.nan),
            "neff": 0.0,
            "count": 0,
        }
    normalized = weight[finite] / total
    selected = values[finite]
    mean = np.sum(normalized[:, None] * selected, axis=0)
    anomaly = selected - mean
    covariance = np.einsum("n,ni,nj->ij", normalized, anomaly, anomaly)
    return {
        "available": True,
        "weight": weight,
        "mean": mean,
        "covariance": covariance,
        "neff": total**2 / max(float(np.sum(weight[finite] ** 2)), 1.0e-30),
        "count": int(np.count_nonzero(finite)),
    }


def _stored_candidate(cube, ti: int, di: int, direction: str):
    """Return one retained-identity-weighted candidate written by Fortran."""
    if direction not in {"up", "down", "combined"}:
        raise ValueError(f"Unsupported PDF-9 candidate direction: {direction}")
    suffix = direction
    required = tuple(
        f"candidate_{name}_{suffix}"
        for name in (
            "valid",
            "w",
            "rt",
            "thl",
            "var_w",
            "var_rt",
            "var_thl",
            "covar_w_rt",
            "covar_w_thl",
            "covar_rt_thl",
        )
    )
    if any(name not in cube for name in required):
        return {
            "available": False,
            "mean": np.full(3, np.nan),
            "covariance": np.full((3, 3), np.nan),
        }
    mean = np.array(
        (
            cube[f"candidate_w_{suffix}"][ti, di],
            cube[f"candidate_rt_{suffix}"][ti, di],
            cube[f"candidate_thl_{suffix}"][ti, di],
        ),
        float,
    )
    covariance = np.array(
        (
            (
                cube[f"candidate_var_w_{suffix}"][ti, di],
                cube[f"candidate_covar_w_rt_{suffix}"][ti, di],
                cube[f"candidate_covar_w_thl_{suffix}"][ti, di],
            ),
            (
                cube[f"candidate_covar_w_rt_{suffix}"][ti, di],
                cube[f"candidate_var_rt_{suffix}"][ti, di],
                cube[f"candidate_covar_rt_thl_{suffix}"][ti, di],
            ),
            (
                cube[f"candidate_covar_w_thl_{suffix}"][ti, di],
                cube[f"candidate_covar_rt_thl_{suffix}"][ti, di],
                cube[f"candidate_var_thl_{suffix}"][ti, di],
            ),
        ),
        float,
    )
    available = bool(
        cube[f"candidate_valid_{suffix}"][ti, di] > 0.5
        and np.all(np.isfinite(mean))
        and np.all(np.isfinite(covariance))
    )
    return {"available": available, "mean": mean, "covariance": covariance}


def _recursive_candidate(cube, ti: int, di: int):
    """Return the pooled bidirectional candidate retained for the next call."""
    return _stored_candidate(cube, ti, di, "combined")


def _active_g3(cube, ti: int, di: int):
    """Return the component-3 geometry and weight used by this PDF call."""
    required = (
        "g3_weight", "g3_w", "g3_rt", "g3_thl", "g3_var_w", "g3_var_rt",
        "g3_var_thl", "g3_corr_w_rt", "g3_corr_w_thl", "g3_corr_rt_thl",
    )
    if not all(name in cube for name in required):
        return {"available": False, "weight": 0.0}
    weight = float(cube["g3_weight"][ti, di])
    mean = np.array(
        (cube["g3_w"][ti, di], cube["g3_rt"][ti, di], cube["g3_thl"][ti, di]),
        float,
    )
    variance = np.array(
        (cube["g3_var_w"][ti, di], cube["g3_var_rt"][ti, di], cube["g3_var_thl"][ti, di]),
        float,
    )
    covariance = np.diag(variance)
    for a, b, name in (
        (0, 1, "g3_corr_w_rt"),
        (0, 2, "g3_corr_w_thl"),
        (1, 2, "g3_corr_rt_thl"),
    ):
        covariance[a, b] = covariance[b, a] = (
            float(cube[name][ti, di])
            * math.sqrt(max(variance[a] * variance[b], 0.0))
        )
    available = bool(
        weight > 1.0e-10
        and np.all(np.isfinite(mean))
        and np.all(np.isfinite(covariance))
    )
    return {
        "available": available,
        "weight": weight if np.isfinite(weight) else 0.0,
        "mean": mean,
        "covariance": covariance,
    }


def _active_residual_component(cube, ti: int, di: int, component: int):
    """Return one actual residual-ADG1 component used by PDF-9."""
    if component not in (1, 2):
        raise ValueError("Residual component must be 1 or 2")
    required = (
        "g12_mixt_frac", "g3_weight", f"g{component}_w", f"g{component}_rt",
        f"g{component}_var_w", f"g{component}_var_rt",
        f"g{component}_corr_w_rt",
    )
    if not all(name in cube for name in required):
        return {"available": False, "weight": 0.0}
    residual_fraction = float(cube["g12_mixt_frac"][ti, di])
    if component == 2:
        residual_fraction = 1.0 - residual_fraction
    g3_weight = _finite_float(cube["g3_weight"][ti, di])
    weight = (1.0 - g3_weight) * residual_fraction
    mean = np.array(
        (cube[f"g{component}_w"][ti, di], cube[f"g{component}_rt"][ti, di]),
        float,
    )
    variance = np.array(
        (cube[f"g{component}_var_w"][ti, di], cube[f"g{component}_var_rt"][ti, di]),
        float,
    )
    correlation = float(cube[f"g{component}_corr_w_rt"][ti, di])
    covariance = np.diag(variance)
    covariance[0, 1] = covariance[1, 0] = (
        np.clip(correlation, -1.0, 1.0)
        * math.sqrt(max(variance[0] * variance[1], 0.0))
    )
    available = bool(
        weight > 1.0e-10
        and np.all(np.isfinite(mean))
        and np.all(np.isfinite(covariance))
        and np.all(variance >= 0.0)
    )
    return {
        "available": available,
        "weight": weight if np.isfinite(weight) else 0.0,
        "mean": mean,
        "covariance": covariance,
    }


def _active_scalar_component(cube, ti: int, di: int, component: int, scalar: str):
    """Return one active component projected onto ``w`` and χ or r_c."""

    if component not in (1, 2, 3) or scalar not in ("chi", "rc"):
        raise ValueError("Scalar component must select G1, G2, or G3 and χ or r_c")
    base = (
        _active_g3(cube, ti, di)
        if component == 3
        else _active_residual_component(cube, ti, di, component)
    )
    fields = (f"g{component}_w", f"g{component}_{scalar}")
    if not base["available"] or not all(field in cube for field in fields):
        return {"available": False, "weight": 0.0}
    mean = np.array((cube[fields[0]][ti, di], cube[fields[1]][ti, di]), float)
    if scalar == "rc":
        return {
            "available": bool(np.all(np.isfinite(mean))),
            "weight": base["weight"],
            "mean": mean,
            "covariance": None,
        }
    variance_fields = (f"g{component}_var_w", f"g{component}_stdev_chi")
    correlation_field = f"g{component}_corr_w_chi"
    if not all(field in cube for field in (*variance_fields, correlation_field)):
        return {"available": False, "weight": 0.0}
    variance = np.array(
        (
            cube[variance_fields[0]][ti, di],
            cube[variance_fields[1]][ti, di] ** 2,
        ),
        float,
    )
    correlation = _finite_float(cube[correlation_field][ti, di])
    covariance = np.diag(variance)
    covariance[0, 1] = covariance[1, 0] = (
        np.clip(correlation, -1.0, 1.0)
        * math.sqrt(max(variance[0] * variance[1], 0.0))
    )
    return {
        "available": bool(
            np.all(np.isfinite(mean))
            and np.all(np.isfinite(covariance))
            and np.all(variance >= 0.0)
        ),
        "weight": base["weight"],
        "mean": mean,
        "covariance": covariance,
    }


def _ellipse(mean, covariance, x_index=1, y_index=0, nrms=2.0, points=121):
    indices = [x_index, y_index]
    matrix = np.asarray(covariance)[np.ix_(indices, indices)]
    if not np.all(np.isfinite(matrix)):
        return np.array([]), np.array([])
    values, vectors = np.linalg.eigh(0.5 * (matrix + matrix.T))
    values = np.maximum(values, 0.0)
    angle = np.linspace(0.0, 2.0 * np.pi, points)
    circle = np.vstack((np.cos(angle), np.sin(angle)))
    path = np.asarray(mean)[indices, None] + vectors @ np.diag(
        nrms * np.sqrt(values)
    ) @ circle
    return path[0], path[1]


def _profile_ensembles(cube, ti: int, scheme: str):
    means = np.full((len(cube["height_m"]), 3), np.nan)
    covariances = np.full((len(cube["height_m"]), 3, 3), np.nan)
    neff = np.zeros(len(cube["height_m"]))
    for di in range(len(cube["height_m"])):
        result = _ensemble(cube, ti, di, scheme)
        if result["available"]:
            means[di] = result["mean"]
            covariances[di] = result["covariance"]
            neff[di] = result["neff"]
    return means, covariances, neff


@lru_cache(maxsize=12)
def _continuity_fields(case_name: str, scheme: str):
    """Vectorize the full time-height reduction once per case and weighting."""
    _meta, cube = _case(case_name)
    valid = cube["parcel_status"] > 0
    parcel_w = np.asarray(cube["parcel_w"], np.float32)
    parcel_rt = np.asarray(cube["parcel_rt"], np.float32)
    if scheme == "uniform":
        weight = np.ones_like(parcel_w)
    elif scheme == "kinetic":
        weight = np.maximum(parcel_w, 0.0) ** 2
    elif scheme == "moist_flux":
        weight = np.maximum(parcel_w, 0.0) * np.maximum(
            parcel_rt - cube["mean_rt"][..., None], 0.0
        )
    elif scheme == "buoyancy_flux":
        weight = np.maximum(parcel_w, 0.0) * np.maximum(
            cube["parcel_buoyancy"], 0.0
        )
    else:
        weight = np.maximum(parcel_w, 0.0)
    weight = np.where(valid & np.isfinite(weight), weight, 0.0)
    total = np.sum(weight, axis=2)
    available = total > 1.0e-24
    safe_total = np.where(available, total, 1.0)
    mean_w = np.sum(weight * parcel_w, axis=2) / safe_total
    mean_rt = np.sum(weight * parcel_rt, axis=2) / safe_total
    neff = total**2 / np.maximum(np.sum(weight**2, axis=2), 1.0e-30)
    mean_w[~available] = np.nan
    mean_rt[~available] = np.nan
    neff[~available] = 0.0
    return mean_w, mean_rt, neff


def _sam_reference(cube, ti: int, reference: str = "all_cloud"):
    index = REFERENCE_NAMES.index(reference)
    return (
        cube["sam_reference_mean"][ti, :, index].astype(float),
        cube["sam_reference_covariance"][ti, :, index].astype(float),
        cube["sam_reference_available"][ti, :, index] > 0,
    )


def _heatmap(snapshot, variable: str = "rt"):
    if variable == "chi":
        scalar = np.asarray(snapshot.chi_samples, float) * 1000.0
        symbol = "χ"
    elif variable == "rc":
        scalar = np.asarray(snapshot.rc_samples, float) * 1000.0
        symbol = "r_c"
    else:
        scalar = np.asarray(snapshot.samples[:, 1], float) * 1000.0
        symbol = "rₜ"
    w = snapshot.samples[:, 0]
    w_prime = w - np.mean(w)
    signed = w_prime * snapshot.rc_samples * 1000.0
    counts, xedges, yedges = np.histogram2d(scalar, w, bins=54)
    transport, _, _ = np.histogram2d(
        scalar, w, bins=(xedges, yedges), weights=signed
    )
    codes, _, _ = signed_bivariate_codes(counts.T, transport.T)
    return go.Heatmap(
        x=0.5 * (xedges[:-1] + xedges[1:]),
        y=0.5 * (yedges[:-1] + yedges[1:]),
        z=codes,
        colorscale=signed_bivariate_reference_colorscale(),
        zmin=0,
        zmax=BIVARIATE_LEVELS * (2 * BIVARIATE_LEVELS - 1) - 1,
        showscale=False,
        hovertemplate=(
            f"{symbol}=%{{x:.3f}} g/kg<br>w=%{{y:.3f}} m/s"
            "<extra>raw SAM bin</extra>"
        ),
    )


def _local_footprint_covariance(cube, ti: int, di: int, scale: float = 0.18):
    """Return a small illustrative w-r_t footprint for parcel centers.

    The ledger contains parcel centers but no within-parcel covariance.  We
    therefore borrow the local CLUBB correlation and a small fraction of its
    width.  This is deliberately distinct from the scientifically scored
    covariance of the full crossing ensemble.
    """
    var_w = max(float(cube["local_var_w"][ti, di]), 0.03**2)
    var_rt = max(float(cube["local_var_rt"][ti, di]), (2.0e-5) ** 2)
    covariance = float(cube["local_cov_w_rt"][ti, di])
    limit = 0.95 * math.sqrt(var_w * var_rt)
    covariance = float(np.clip(covariance, -limit, limit))
    result = np.zeros((3, 3), float)
    result[0, 0] = var_w * scale**2
    result[1, 1] = var_rt * scale**2
    result[0, 1] = result[1, 0] = covariance * scale**2
    return result


def _add_directional_parcel_footprints(
    figure,
    cube,
    ti: int,
    di: int,
    result,
    *,
    direction: str,
):
    downward = direction == "down"
    prefix = "down_" if downward else ""
    color = COLORS["down"] if downward else COLORS["up"]
    rgba = "rgba(251,113,133,0.48)" if downward else "rgba(34,211,238,0.48)"
    label = "Downward" if downward else "Individual"
    valid = (
        (cube[f"{prefix}parcel_status"][ti, di] > 0)
        & np.isfinite(cube[f"{prefix}parcel_w"][ti, di])
        & np.isfinite(cube[f"{prefix}parcel_rt"][ti, di])
    )
    launches = np.where(valid)[0]
    if not launches.size:
        return
    covariance = _local_footprint_covariance(cube, ti, di)
    x_path: list[float | None] = []
    y_path: list[float | None] = []
    for launch in launches:
        mean = np.array(
            (
                cube[f"{prefix}parcel_w"][ti, di, launch],
                cube[f"{prefix}parcel_rt"][ti, di, launch],
                cube[f"{prefix}parcel_thl"][ti, di, launch],
            ),
            float,
        )
        x, y = _ellipse(mean, covariance, points=45)
        x_path.extend((x * 1000.0).tolist() + [None])
        y_path.extend(y.tolist() + [None])
    figure.add_trace(
        go.Scattergl(
            x=x_path,
            y=y_path,
            mode="lines",
            line={"color": rgba, "width": 1},
            name=f"{label} parcel footprints (illustrative local width)",
            hoverinfo="skip",
        )
    )
    raw_weight = np.asarray(result["weight"], float)[launches]
    maximum = max(float(np.nanmax(raw_weight)), 1.0e-30)
    figure.add_trace(
        go.Scattergl(
            x=cube[f"{prefix}parcel_rt"][ti, di, launches] * 1000.0,
            y=cube[f"{prefix}parcel_w"][ti, di, launches],
            mode="markers",
            marker={
                "color": color,
                "size": 4.0 + 5.0 * np.sqrt(np.maximum(raw_weight, 0.0) / maximum),
                "opacity": 0.78,
                "symbol": "triangle-down" if downward else "circle",
            },
            customdata=np.column_stack(
                (cube["launch_height_m"][launches], raw_weight)
            ),
            hovertemplate=(
                "launch=%{customdata[0]:.0f} m<br>"
                "r_t=%{x:.3f} g/kg<br>w=%{y:.3f} m/s<br>"
                "weight=%{customdata[1]:.3g}<extra>parcel crossing</extra>"
            ),
            name=f"{label} parcel crossing centers",
        )
    )


def _add_parcel_footprints(figure, cube, ti: int, di: int, result):
    _add_directional_parcel_footprints(
        figure, cube, ti, di, result, direction="up"
    )


def _add_downward_parcel_footprints(figure, cube, ti: int, di: int, result):
    _add_directional_parcel_footprints(
        figure, cube, ti, di, result, direction="down"
    )


def _add_occupancy_center(figure, cube, ti: int, di: int):
    """Overlay the equal-crossing center independently of the selected ellipse."""
    occupancy = _ensemble(cube, ti, di, "uniform")
    if occupancy["available"]:
        figure.add_trace(
            go.Scatter(
                x=[occupancy["mean"][1] * 1000.0],
                y=[occupancy["mean"][0]],
                mode="markers",
                marker={
                    "color": COLORS["purple"],
                    "size": 15,
                    "symbol": "diamond-open",
                    "line": {"color": COLORS["purple"], "width": 3},
                },
                name="Occupancy-weighted candidate center",
                hovertemplate=(
                    "r_t=%{x:.3f} g/kg<br>w=%{y:.3f} m/s"
                    "<extra>equal vote per crossing</extra>"
                ),
            )
        )
    return occupancy


def _scalar_overlay_figure(cube, ti, di, snapshot, scalar: str):
    """Draw the raw-SAM scalar heatmap and its active PDF-component geometry."""

    label = OVERLAY_VARIABLES[scalar]["label"]
    figure = go.Figure(_heatmap(snapshot, scalar))
    for component in (1, 2, 3):
        state = _active_scalar_component(cube, ti, di, component, scalar)
        if not state["available"]:
            continue
        color = COLORS["recursive"] if component == 3 else COLORS[f"g{component}"]
        marker = "star" if component == 3 else "x"
        if state["covariance"] is not None:
            x, y = _ellipse(state["mean"], state["covariance"])
            figure.add_trace(
                go.Scatter(
                    x=x * 1000.0,
                    y=y,
                    mode="lines",
                    line={"color": color, "width": 5 if component == 3 else 4, "dash": "solid" if component == 3 else "dash"},
                    name=(
                        f"Active Fortran G{component} · 2σ "
                        f"(a{component}={state['weight']:.3f})"
                    ),
                )
            )
        figure.add_trace(
            go.Scatter(
                x=[state["mean"][1] * 1000.0],
                y=[state["mean"][0]],
                mode="markers",
                marker={"color": color, "size": 17 if component == 3 else 14, "symbol": marker},
                name=(
                    f"Active Fortran G{component} center"
                    if scalar == "chi"
                    else f"Active Fortran G{component} r_c center"
                ),
                hovertemplate=(
                    f"{label}=%{{x:.3f}} g/kg<br>w=%{{y:.3f}} m/s"
                    f"<extra>G{component}</extra>"
                ),
            )
        )
    caveat = (
        "PDF component geometry over raw SAM"
        if scalar == "chi"
        else "Raw SAM cloud-water structure with PDF component centers"
    )
    figure.update_layout(
        template="plotly_dark",
        title=caveat,
        xaxis_title=OVERLAY_VARIABLES[scalar]["axis"],
        yaxis_title="Vertical velocity w [m/s]",
        height=610,
        margin={"l": 60, "r": 20, "t": 65, "b": 60},
        legend={"orientation": "h", "y": -0.18},
        uirevision=f"pdf9-overlay-{scalar}",
    )
    if scalar == "rc":
        figure.add_annotation(
            text="r_c is nonlinear; component-level w–r_c covariances are not emitted, so centers are shown without Gaussian contours.",
            xref="paper",
            yref="paper",
            x=0.5,
            y=-0.29,
            showarrow=False,
            font={"color": "#94a3b8", "size": 11},
        )
    return figure


def overlay_figure(cube, meta, ti, di, scheme, variable: str = "rt"):
    result = _ensemble(cube, ti, di, scheme)
    downward_weights = _downward_ensemble(cube, ti, di)
    downward = _stored_candidate(cube, ti, di, "down")
    recursive = _recursive_candidate(cube, ti, di)
    residual_components = (
        _active_residual_component(cube, ti, di, 1),
        _active_residual_component(cube, ti, di, 2),
    )
    active_g3 = _active_g3(cube, ti, di)
    seconds = int(round(float(cube["sam_time_seconds"][ti])))
    height = float(cube["sam_height_m"][di])
    snapshot = _raw_snapshot(meta["sam_run_dir"], seconds, height)
    if variable in ("chi", "rc"):
        return _scalar_overlay_figure(cube, ti, di, snapshot, variable), snapshot, result
    figure = go.Figure(_heatmap(snapshot, "rt"))
    _add_parcel_footprints(figure, cube, ti, di, result)
    _add_downward_parcel_footprints(
        figure, cube, ti, di, downward_weights
    )
    if result["available"]:
        x, y = _ellipse(result["mean"], result["covariance"])
        figure.add_trace(
            go.Scatter(x=x * 1000.0, y=y, mode="lines", line={"color": COLORS["candidate"], "width": 3}, name="Parcel crossing ensemble · 2σ")
        )
        if scheme != "uniform":
            figure.add_trace(
                go.Scatter(x=[result["mean"][1] * 1000.0], y=[result["mean"][0]], mode="markers", marker={"color": COLORS["candidate"], "size": 13, "symbol": "cross"}, name="Selected-weighting candidate center")
            )
    _add_occupancy_center(figure, cube, ti, di)
    for component, residual in enumerate(residual_components, start=1):
        if not residual["available"]:
            continue
        color = COLORS[f"g{component}"]
        x, y = _ellipse(residual["mean"], residual["covariance"])
        figure.add_trace(
            go.Scatter(
                x=x * 1000.0,
                y=y,
                mode="lines",
                line={"color": color, "width": 4, "dash": "dash"},
                name=(
                    f"Active Fortran G{component} · 2σ "
                    f"(a{component}={residual['weight']:.3f})"
                ),
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[residual["mean"][1] * 1000.0],
                y=[residual["mean"][0]],
                mode="markers",
                marker={"color": color, "size": 14, "symbol": "x"},
                name=f"Active Fortran G{component} center",
            )
        )
    if active_g3["available"]:
        x, y = _ellipse(active_g3["mean"], active_g3["covariance"])
        figure.add_trace(
            go.Scatter(
                x=x * 1000.0,
                y=y,
                mode="lines",
                line={"color": COLORS["recursive"], "width": 5},
                name=f"Active Fortran G3 · 2σ (a₃={active_g3['weight']:.3f})",
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[active_g3["mean"][1] * 1000.0],
                y=[active_g3["mean"][0]],
                mode="markers",
                marker={
                    "color": COLORS["recursive"],
                    "size": 17,
                    "symbol": "star",
                },
                name="Active Fortran G3 center",
            )
        )
    if recursive["available"]:
        x, y = _ellipse(recursive["mean"], recursive["covariance"])
        figure.add_trace(
            go.Scatter(
                x=x * 1000.0,
                y=y,
                mode="lines",
                line={"color": COLORS["recursive"], "width": 4},
                name="Next pooled G3 candidate · 2σ",
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[recursive["mean"][1] * 1000.0],
                y=[recursive["mean"][0]],
                mode="markers",
                marker={
                    "color": COLORS["recursive"],
                    "size": 17,
                    "symbol": "star-open",
                    "line": {"color": COLORS["recursive"], "width": 3},
                },
                name="Next pooled G3 center",
            )
        )
    if downward["available"]:
        x, y = _ellipse(downward["mean"], downward["covariance"])
        figure.add_trace(
            go.Scatter(
                x=x * 1000.0,
                y=y,
                mode="lines",
                line={"color": "#fb7185", "width": 3, "dash": "dash"},
                name="Fortran downward candidate · 2σ",
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[downward["mean"][1] * 1000.0],
                y=[downward["mean"][0]],
                mode="markers",
                marker={"color": "#fb7185", "size": 13, "symbol": "triangle-down"},
                name="Fortran downward candidate center",
            )
        )
    sam_mean, sam_cov, available = _sam_reference(cube, ti)
    if available[di]:
        x, y = _ellipse(sam_mean[di], sam_cov[di])
        figure.add_trace(go.Scatter(x=x * 1000.0, y=y, mode="lines", line={"color": COLORS["sam"], "width": 2, "dash": "dash"}, name="SAM rᶜ-weighted shape · 2σ"))
        figure.add_trace(go.Scatter(x=[sam_mean[di, 1] * 1000.0], y=[sam_mean[di, 0]], mode="markers", marker={"color": COLORS["sam"], "size": 12, "symbol": "star"}, name="SAM rᶜ-weighted center"))
    figure.update_layout(
        template="plotly_dark",
        title=f"Candidate transport Gaussian over raw SAM · {seconds/60:.0f} min / {height:.0f} m",
        xaxis_title="Total water rₜ [g/kg]",
        yaxis_title="Vertical velocity w [m/s]",
        height=610,
        margin={"l": 60, "r": 20, "t": 65, "b": 60},
        legend={"orientation": "h", "y": -0.18},
        uirevision=f"pdf9-overlay-{meta['name']}",
    )
    return figure, snapshot, result


def component_weights_figure(cube, ti: int):
    """Show the three actual PDF mixture weights through the selected column."""

    required = ("g12_mixt_frac", "g3_weight", "height_m")
    if not all(name in cube for name in required):
        return _empty("PDF component weights", "Regenerate the PDF-9 cache")
    g3 = np.clip(np.asarray(cube["g3_weight"][ti], float), 0.0, 1.0)
    g12 = np.clip(np.asarray(cube["g12_mixt_frac"][ti], float), 0.0, 1.0)
    residual = 1.0 - g3
    weights = (
        ("G1", residual * g12, COLORS["g1"]),
        ("G2", residual * (1.0 - g12), COLORS["g2"]),
        ("G3", g3, COLORS["recursive"]),
    )
    figure = go.Figure()
    for name, values, color in weights:
        figure.add_trace(
            go.Scatter(
                x=values,
                y=cube["height_m"],
                mode="lines+markers",
                line={"color": color, "width": 3},
                marker={"size": 5},
                name=name,
                hovertemplate=(
                    f"{name} weight=%{{x:.4f}}<br>z=%{{y:.0f}} m"
                    f"<extra>{name}</extra>"
                ),
                connectgaps=False,
            )
        )
    figure.update_layout(
        template="plotly_dark",
        height=530,
        title=(
            "Actual PDF component weights through the column · "
            f"{cube['time_seconds'][ti] / 60.0:.0f} min"
        ),
        xaxis={"title": "PDF mixture weight", "range": [0.0, 1.0]},
        yaxis={"title": "Altitude [m]"},
        legend={"orientation": "h", "y": -0.18},
        uirevision="pdf9-component-weights",
    )
    return figure


def tracks_figure(cube, ti, di):
    figure = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=(
            "Parcel total-water tracks",
            "Parcel liquid-water potential-temperature tracks",
            "Parcel vertical-speed tracks",
        ),
    )
    status = cube["parcel_status"][ti] > 0
    candidates = np.where(np.any(status, axis=0))[0]
    if candidates.size:
        rt_path: list[float | None] = []
        thl_path: list[float | None] = []
        w_path: list[float | None] = []
        z_path: list[float | None] = []
        for launch in candidates:
            valid = status[:, launch]
            z = cube["height_m"][valid]
            if not len(z):
                continue
            # The first stored ledger point is the first crossed interface.
            # Prepending the launch height at that same state makes the origin
            # of every trajectory explicit without inventing a launch state.
            z = np.concatenate(([cube["launch_height_m"][launch]], z))
            rt = np.asarray(cube["parcel_rt"][ti, valid, launch], float) * 1000.0
            thl = np.asarray(cube["parcel_thl"][ti, valid, launch], float)
            w = np.asarray(cube["parcel_w"][ti, valid, launch], float)
            rt = np.concatenate(([rt[0]], rt))
            thl = np.concatenate(([thl[0]], thl))
            w = np.concatenate(([w[0]], w))
            rt_path.extend(rt.tolist() + [None])
            thl_path.extend(thl.tolist() + [None])
            w_path.extend(w.tolist() + [None])
            z_path.extend(z.tolist() + [None])
        line = {"color": "rgba(34,211,238,0.34)", "width": 1.2}
        figure.add_trace(
            go.Scattergl(
                x=rt_path,
                y=z_path,
                mode="lines",
                name=f"All {len(candidates)} upward launch parcels",
                line=line,
                hoverinfo="skip",
            ),
            row=1,
            col=1,
        )
        figure.add_trace(
            go.Scattergl(
                x=thl_path,
                y=z_path,
                mode="lines",
                name="All upward launch parcels",
                line=line,
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
        figure.add_trace(
            go.Scattergl(
                x=w_path,
                y=z_path,
                mode="lines",
                name="All upward launch parcels",
                line=line,
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1,
            col=3,
        )

        selected = np.where(status[di])[0]
        if selected.size:
            figure.add_trace(
                go.Scattergl(
                    x=cube["parcel_rt"][ti, di, selected] * 1000.0,
                    y=np.full(selected.size, cube["height_m"][di]),
                    mode="markers",
                    marker={"color": COLORS["candidate"], "size": 7},
                    name="Crossings at selected altitude",
                ),
                row=1,
                col=1,
            )
            figure.add_trace(
                go.Scattergl(
                    x=cube["parcel_thl"][ti, di, selected],
                    y=np.full(selected.size, cube["height_m"][di]),
                    mode="markers",
                    marker={"color": COLORS["candidate"], "size": 7},
                    showlegend=False,
                    name="Crossings at selected altitude",
                ),
                row=1,
                col=2,
            )
            figure.add_trace(
                go.Scattergl(
                    x=cube["parcel_w"][ti, di, selected],
                    y=np.full(selected.size, cube["height_m"][di]),
                    mode="markers",
                    marker={"color": COLORS["candidate"], "size": 7},
                    showlegend=False,
                    name="Crossings at selected altitude",
                ),
                row=1,
                col=3,
            )

    down_status = cube["down_parcel_status"][ti] > 0
    down_candidates = np.where(np.any(down_status, axis=0))[0]
    if down_candidates.size:
        rt_path: list[float | None] = []
        thl_path: list[float | None] = []
        w_path: list[float | None] = []
        z_path: list[float | None] = []
        for launch in down_candidates:
            # Destination storage is bottom-to-top for both ledgers.  A
            # sinking trajectory must instead be drawn from its launch to the
            # nearest lower crossing and then continue downward.  Reversing
            # the valid destination indices prevents the apparent "bounce"
            # back upward seen when the raw storage order is connected.
            destinations = np.where(down_status[:, launch])[0][::-1]
            if not destinations.size:
                continue
            z = np.asarray(cube["height_m"][destinations], float)
            z = np.concatenate(([cube["launch_height_m"][launch]], z))
            rt = np.asarray(
                cube["down_parcel_rt"][ti, destinations, launch], float
            ) * 1000.0
            thl = np.asarray(
                cube["down_parcel_thl"][ti, destinations, launch], float
            )
            w = np.asarray(
                cube["down_parcel_w"][ti, destinations, launch], float
            )
            rt_path.extend(np.concatenate(([rt[0]], rt)).tolist() + [None])
            thl_path.extend(np.concatenate(([thl[0]], thl)).tolist() + [None])
            w_path.extend(np.concatenate(([w[0]], w)).tolist() + [None])
            z_path.extend(z.tolist() + [None])
        line = {"color": "rgba(251,113,133,0.38)", "width": 1.2}
        for column, path, name in (
            (1, rt_path, f"All {len(down_candidates)} downward launch parcels"),
            (2, thl_path, "All downward launch parcels"),
            (3, w_path, "All downward launch parcels"),
        ):
            figure.add_trace(
                go.Scattergl(
                    x=path,
                    y=z_path,
                    mode="lines",
                    name=name,
                    line=line,
                    hoverinfo="skip",
                    showlegend=column == 1,
                ),
                row=1,
                col=column,
            )

        selected = np.where(down_status[di])[0]
        if selected.size:
            for column, field, scale in (
                (1, "down_parcel_rt", 1000.0),
                (2, "down_parcel_thl", 1.0),
                (3, "down_parcel_w", 1.0),
            ):
                figure.add_trace(
                    go.Scattergl(
                        x=cube[field][ti, di, selected] * scale,
                        y=np.full(selected.size, cube["height_m"][di]),
                        mode="markers",
                        marker={"color": "#fb7185", "size": 7, "symbol": "triangle-down"},
                        name="Downward crossings at selected altitude",
                        showlegend=column == 1,
                    ),
                    row=1,
                    col=column,
                )

    mean_rt = np.asarray(cube["mean_rt"][ti], float) * 1000.0
    sigma_rt = np.sqrt(np.maximum(cube["local_var_rt"][ti], 0.0)) * 1000.0
    z = np.asarray(cube["height_m"], float)
    figure.add_trace(
        go.Scatter(
            x=np.concatenate((mean_rt + 2.0 * sigma_rt, (mean_rt - 2.0 * sigma_rt)[::-1])),
            y=np.concatenate((z, z[::-1])),
            mode="lines",
            line={"width": 0},
            fill="toself",
            fillcolor="rgba(148,163,184,0.16)",
            name="CLUBB grid mean r_t +/- 2 sigma",
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )
    sam, _, available = _sam_reference(cube, ti)
    clear_moist, _, clear_moist_available = _sam_reference(
        cube, ti, "clear_moist_up"
    )
    clear_moist = np.where(clear_moist_available[:, None], clear_moist, np.nan)
    figure.add_trace(go.Scatter(x=mean_rt, y=cube["height_m"], mode="lines", line={"color": COLORS["mean"], "dash": "dash", "width": 2}, name="CLUBB grid mean"), row=1, col=1)
    figure.add_trace(go.Scatter(x=cube["mean_thl"][ti], y=cube["height_m"], mode="lines", line={"color": COLORS["mean"], "dash": "dash", "width": 2}, name="CLUBB grid mean", showlegend=False), row=1, col=2)
    figure.add_trace(go.Scatter(x=sam[:, 1] * 1000.0, y=cube["height_m"], mode="lines", line={"color": COLORS["sam"], "width": 3}, name="SAM cloud-water-weighted center", connectgaps=False), row=1, col=1)
    figure.add_trace(go.Scatter(x=sam[:, 2], y=cube["height_m"], mode="lines", line={"color": COLORS["sam"], "width": 3}, name="SAM cloud-water-weighted center", showlegend=False, connectgaps=False), row=1, col=2)
    figure.add_trace(go.Scatter(x=sam[:, 0], y=cube["height_m"], mode="lines", line={"color": COLORS["sam"], "width": 3}, name="SAM cloud-water-weighted center", showlegend=False, connectgaps=False), row=1, col=3)
    for column, state_index, scale in ((1, 1, 1000.0), (2, 2, 1.0), (3, 0, 1.0)):
        figure.add_trace(
            go.Scatter(
                x=clear_moist[:, state_index] * scale,
                y=cube["height_m"],
                mode="lines",
                line={"color": COLORS["up"], "width": 3, "dash": "dot"},
                name="SAM clear moist updraft",
                showlegend=column == 1,
                connectgaps=False,
            ),
            row=1,
            col=column,
        )
    if "candidate_valid_combined" in cube:
        valid_candidate = (
            np.asarray(cube["candidate_valid_combined"][ti], float) > 0.5
        )
        candidate_fields = (
            (1, "candidate_rt_combined", 1000.0),
            (2, "candidate_thl_combined", 1.0),
            (3, "candidate_w_combined", 1.0),
        )
        for column, field, scale in candidate_fields:
            values = np.where(
                valid_candidate,
                np.asarray(cube[field][ti], float) * scale,
                np.nan,
            )
            figure.add_trace(
                go.Scatter(
                    x=values,
                    y=cube["height_m"],
                    mode="lines+markers",
                    line={"color": COLORS["recursive"], "width": 3},
                    marker={"size": 5, "symbol": "star-open"},
                    name="Fortran pooled G3 candidate",
                    showlegend=column == 1,
                    connectgaps=False,
                ),
                row=1,
                col=column,
            )
    figure.add_hline(y=float(cube["height_m"][di]), line_color="#64748b", line_dash="dot")
    figure.update_xaxes(title_text="rₜ [g/kg]", row=1, col=1)
    figure.update_xaxes(title_text="θₗ [K]", row=1, col=2)
    figure.update_xaxes(title_text="w [m/s]", row=1, col=3)
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(template="plotly_dark", height=610, title="Every upward and downward parcel trajectory through the column", legend={"orientation": "h", "y": -0.18}, uirevision="pdf9-tracks")
    return figure


def reach_and_support_figure(cube, ti, di, scheme):
    z = cube["height_m"]
    step = max(1, len(z) // 28)
    nz = len(z)
    contribution = np.zeros((nz, len(cube["launch_height_m"])))
    neff = np.zeros(nz)
    for destination in range(nz):
        result = _ensemble(cube, ti, destination, scheme)
        if np.sum(result["weight"]) > 0:
            contribution[destination] = result["weight"] / np.sum(result["weight"])
        neff[destination] = result["neff"]

    figure = make_subplots(
        rows=1,
        cols=3,
        column_widths=(0.38, 0.42, 0.20),
        subplot_titles=(
            "Upward and downward mixing reach",
            "Normalized launch contribution",
            "Effective parcel count",
        ),
        shared_yaxes=False,
    )
    for field, color, sign, name in (
        ("lscale_up", COLORS["up"], 1.0, "Upward reach"),
        ("lscale_down", COLORS["down"], -1.0, "Downward reach"),
    ):
        for index in range(0, len(z), step):
            end = z[index] + sign * float(cube[field][ti, index])
            figure.add_trace(
                go.Scatter(
                    x=[z[index], z[index]],
                    y=[z[index], end],
                    mode="lines+markers",
                    line={"color": color, "width": 2},
                    marker={"size": [5, 8], "symbol": ["square-open", "arrow-up" if sign > 0 else "arrow-down"]},
                    name=name,
                    legendgroup=name,
                    showlegend=index == 0,
                    hovertemplate=(
                        f"launch={z[index]:.0f} m<br>arrival=%{{y:.0f}} m"
                        f"<br>reach={abs(end-z[index]):.0f} m<extra>{name}</extra>"
                    ),
                ),
                row=1,
                col=1,
            )
    figure.add_trace(
        go.Heatmap(
            x=cube["launch_height_m"],
            y=z,
            z=contribution,
            colorscale="Viridis",
            colorbar={"title": "q/Σq", "len": 0.78},
            showscale=True,
        ),
        row=1,
        col=2,
    )
    figure.add_trace(
        go.Scatter(
            x=neff,
            y=z,
            mode="lines",
            line={"color": COLORS["candidate"], "width": 3},
            name="N_eff",
        ),
        row=1,
        col=3,
    )
    for column in (1, 2, 3):
        figure.add_hline(
            y=float(z[di]), line_color=COLORS["sam"], line_dash="dot", row=1, col=column
        )
    figure.update_xaxes(title_text="Launch altitude [m]", row=1, col=1)
    figure.update_xaxes(title_text="Launch altitude [m]", row=1, col=2)
    figure.update_xaxes(title_text="N_eff", row=1, col=3)
    figure.update_yaxes(title_text="Arrival altitude [m]", row=1, col=1)
    figure.update_yaxes(title_text="Destination altitude [m]", row=1, col=2)
    figure.update_yaxes(title_text="Destination altitude [m]", row=1, col=3)
    figure.update_layout(
        template="plotly_dark",
        height=570,
        title="Directional mixing reach and source support",
        legend={"orientation": "h", "y": -0.18},
        uirevision="pdf9-reach-and-support",
    )
    return figure


def entrainment_environment_figure(cube, ti, di):
    """Show the exact scalar environment seen by each directional parcel."""
    z = np.asarray(cube["height_m"], float)
    mean_rt = np.asarray(cube["mean_rt"][ti], float) * 1000.0
    mean_thl = np.asarray(cube["mean_thl"][ti], float)
    g3_weight = np.asarray(cube["g3_weight"][ti], float)
    applied_weight = np.asarray(cube["entrain_g3_weight"][ti], float)
    active = applied_weight > 0.0

    figure = make_subplots(
        rows=1,
        cols=3,
        shared_yaxes=True,
        column_widths=(0.39, 0.39, 0.22),
        subplot_titles=(
            "Total-water entrainment target",
            "Temperature entrainment target",
            "Applied G3 fraction",
        ),
    )
    for column, mean, up_field, down_field, scale, label in (
        (1, mean_rt, "entrain_rt_up", "entrain_rt_down", 1000.0, "rₜ"),
        (2, mean_thl, "entrain_thl_up", "entrain_thl_down", 1.0, "θₗ"),
    ):
        figure.add_trace(
            go.Scatter(
                x=mean,
                y=z,
                mode="lines",
                line={"color": COLORS["mean"], "dash": "dash", "width": 2},
                name="Grid-mean environment",
                showlegend=column == 1,
            ),
            row=1,
            col=column,
        )
        for field, color, name in (
            (up_field, COLORS["up"], "Upward parcel target"),
            (down_field, COLORS["down"], "Downward parcel target"),
        ):
            values = np.asarray(cube[field][ti], float) * scale
            figure.add_trace(
                go.Scatter(
                    x=np.where(active, values, np.nan),
                    y=z,
                    mode="lines+markers",
                    line={"color": color, "width": 3},
                    marker={"size": 5},
                    name=name,
                    legendgroup=name,
                    showlegend=column == 1,
                    connectgaps=False,
                    hovertemplate=(
                        f"z=%{{y:.0f}} m<br>{label}=%{{x:.5g}}<extra>{name}</extra>"
                    ),
                ),
                row=1,
                col=column,
            )
    figure.add_trace(
        go.Scatter(
            x=g3_weight,
            y=z,
            mode="lines",
            line={"color": COLORS["purple"], "dash": "dot", "width": 2},
            name="PDF G3 weight",
        ),
        row=1,
        col=3,
    )
    figure.add_trace(
        go.Scatter(
            x=applied_weight,
            y=z,
            mode="lines+markers",
            line={"color": COLORS["recursive"], "width": 3},
            marker={"size": 5},
            name="Weight used by entrainment",
        ),
        row=1,
        col=3,
    )
    figure.add_hline(y=float(z[di]), line_color=COLORS["sam"], line_width=2)
    figure.update_xaxes(title_text="rₜ [g/kg]", row=1, col=1)
    figure.update_xaxes(title_text="θₗ [K]", row=1, col=2)
    figure.update_xaxes(title_text="Mixture fraction", row=1, col=3)
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(
        template="plotly_dark",
        height=570,
        title=(
            f"G3-aware scalar entrainment environment · "
            f"{cube['time_seconds'][ti] / 60.0:.0f} min"
        ),
        legend={"orientation": "h", "y": -0.18},
        uirevision="pdf9-entrainment-environment",
    )
    return figure


def sources_figure(cube, ti, di, scheme):
    nz = len(cube["height_m"])
    matrix = np.zeros((nz, len(cube["launch_height_m"])))
    neff = np.zeros(nz)
    for destination in range(nz):
        result = _ensemble(cube, ti, destination, scheme)
        weight = result["weight"]
        if np.sum(weight) > 0:
            matrix[destination] = weight / np.sum(weight)
        neff[destination] = result["neff"]
    figure = make_subplots(rows=1, cols=2, column_widths=(0.72, 0.28), subplot_titles=("Normalized launch contribution", "Effective parcel count"))
    figure.add_trace(go.Heatmap(x=cube["launch_height_m"], y=cube["height_m"], z=matrix, colorscale="Viridis", colorbar={"title": "q/Σq", "len": 0.8}), row=1, col=1)
    figure.add_trace(go.Scatter(x=neff, y=cube["height_m"], mode="lines", line={"color": COLORS["candidate"], "width": 3}, name="N_eff"), row=1, col=2)
    figure.add_hline(y=float(cube["height_m"][di]), line_color=COLORS["sam"], line_dash="dot")
    figure.update_xaxes(title_text="Launch altitude [m]", row=1, col=1)
    figure.update_xaxes(title_text="N_eff", row=1, col=2)
    figure.update_yaxes(title_text="Destination altitude [m]", row=1, col=1)
    figure.update_layout(template="plotly_dark", height=600, title="Who built the candidate at each destination?", uirevision="pdf9-source-map")
    return figure


def error_figure(cube, ti, scheme):
    means, covariances, neff = _profile_ensembles(cube, ti, scheme)
    sam_mean, sam_cov, available = _sam_reference(cube, ti)
    z = cube["height_m"]
    sigma = np.sqrt(np.maximum(np.diagonal(sam_cov, axis1=1, axis2=2), 1.0e-30))
    error = (means - sam_mean) / sigma
    candidate_corr = covariances[:, 0, 1] / np.sqrt(np.maximum(covariances[:, 0, 0] * covariances[:, 1, 1], 1.0e-30))
    sam_corr = sam_cov[:, 0, 1] / np.sqrt(np.maximum(sam_cov[:, 0, 0] * sam_cov[:, 1, 1], 1.0e-30))
    figure = make_subplots(rows=1, cols=3, subplot_titles=("Center errors", "w–rₜ correlation", "Support"), shared_yaxes=True)
    for index, label, color in ((0, "w", COLORS["up"]), (1, "rₜ", COLORS["sam"]), (2, "θₗ", COLORS["down"])):
        figure.add_trace(go.Scatter(x=np.where(available, error[:, index], np.nan), y=z, mode="lines", name=f"{label} error [SAM σ]", line={"color": color}), row=1, col=1)
    figure.add_trace(go.Scatter(x=candidate_corr, y=z, mode="lines", name="candidate ρ(w,rₜ)", line={"color": COLORS["candidate"], "width": 3}), row=1, col=2)
    figure.add_trace(go.Scatter(x=sam_corr, y=z, mode="lines", name="SAM cloud-tail ρ(w,rₜ)", line={"color": COLORS["sam"], "dash": "dash"}), row=1, col=2)
    figure.add_trace(go.Scatter(x=neff, y=z, mode="lines", name="N_eff", line={"color": COLORS["purple"], "width": 3}), row=1, col=3)
    figure.update_xaxes(title_text="standard deviations", row=1, col=1)
    figure.update_xaxes(title_text="correlation", row=1, col=2, range=[-1.05, 1.05])
    figure.update_xaxes(title_text="effective launches", row=1, col=3)
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(template="plotly_dark", height=570, title="Center, shape, and sampling error against raw SAM", legend={"orientation": "h", "y": -0.18}, uirevision="pdf9-errors")
    return figure


def anomaly_figure(cube, ti, di):
    figure = make_subplots(rows=1, cols=3, subplot_titles=("Moisture anomaly", "θₗ anomaly", "Vertical speed"), shared_yaxes=True)
    active = np.where(cube["parcel_status"][ti, di] > 0)[0]
    if not active.size:
        return _empty("Parcel anomaly evolution")
    selected = active[np.unique(np.rint(np.linspace(0, len(active) - 1, min(6, len(active)))).astype(int))]
    for launch in selected:
        valid = cube["parcel_status"][ti, :, launch] > 0
        z = cube["height_m"][valid]
        label = f"{cube['launch_height_m'][launch]:.0f} m launch"
        figure.add_trace(go.Scatter(x=(cube["parcel_rt"][ti, valid, launch] - cube["mean_rt"][ti, valid]) * 1000.0, y=z, mode="lines", name=label, legendgroup=str(launch)), row=1, col=1)
        figure.add_trace(go.Scatter(x=cube["parcel_thl"][ti, valid, launch] - cube["mean_thl"][ti, valid], y=z, mode="lines", name=label, legendgroup=str(launch), showlegend=False), row=1, col=2)
        figure.add_trace(go.Scatter(x=cube["parcel_w"][ti, valid, launch], y=z, mode="lines", name=label, legendgroup=str(launch), showlegend=False), row=1, col=3)
    figure.add_vline(x=0, line_color=COLORS["mean"], line_dash="dash")
    figure.update_xaxes(title_text="rₜ parcel − environment [g/kg]", row=1, col=1)
    figure.update_xaxes(title_text="θₗ parcel − environment [K]", row=1, col=2)
    figure.update_xaxes(title_text="w [m/s]", row=1, col=3)
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(template="plotly_dark", height=560, title="How entrainment erodes each surviving parcel anomaly", legend={"orientation": "h", "y": -0.2}, uirevision="pdf9-anomalies")
    return figure


def bakeoff_figure(cube, ti):
    sam_mean, sam_cov, available = _sam_reference(cube, ti)
    sigma = np.sqrt(np.maximum(np.diagonal(sam_cov, axis1=1, axis2=2), 1.0e-30))
    labels, total, w_error, rt_error, shape_error = [], [], [], [], []
    for scheme, label in WEIGHTINGS.items():
        means, covariance, neff = _profile_ensembles(cube, ti, scheme)
        valid = available & (neff >= 1.0)
        center = (means - sam_mean) / sigma
        wrt = covariance[:, 0, 1] / np.sqrt(np.maximum(covariance[:, 0, 0] * covariance[:, 1, 1], 1.0e-30))
        sam_wrt = sam_cov[:, 0, 1] / np.sqrt(np.maximum(sam_cov[:, 0, 0] * sam_cov[:, 1, 1], 1.0e-30))
        def rms(values):
            selected = np.asarray(values)[valid]
            return float(np.sqrt(np.nanmean(selected**2))) if np.any(valid) else math.nan
        labels.append(label)
        w_error.append(rms(center[:, 0]))
        rt_error.append(rms(center[:, 1]))
        shape_error.append(rms(wrt - sam_wrt))
        total.append(np.nanmean((center[valid, 0] ** 2 + center[valid, 1] ** 2 + (wrt[valid] - sam_wrt[valid]) ** 2) / 3.0) ** 0.5 if np.any(valid) else math.nan)
    figure = go.Figure()
    for values, name, color in ((w_error, "w center [σ]", COLORS["up"]), (rt_error, "rₜ center [σ]", COLORS["sam"]), (shape_error, "ρ(w,rₜ) error", COLORS["purple"]), (total, "combined", COLORS["candidate"])):
        figure.add_bar(x=labels, y=values, name=name, marker_color=color)
    figure.update_layout(template="plotly_dark", barmode="group", height=510, title="Weighting comparison over this entire SAM column", yaxis_title="lower is better", legend={"orientation": "h", "y": -0.28}, uirevision="pdf9-bakeoff")
    return figure


@lru_cache(maxsize=16)
def _whole_case_center_metrics(case_name: str, scheme: str):
    """Score a dashboard-only crossing reduction over every cached SAM plane."""
    _meta, cube = _case(case_name)
    candidate_w, candidate_rt, neff = _continuity_fields(case_name, scheme)
    sam_mean, _sam_covariance, sam_available = _sam_reference(cube, slice(None))
    valid = (
        sam_available
        & (neff >= 1.0)
        & np.isfinite(candidate_w)
        & np.isfinite(candidate_rt)
        & np.isfinite(sam_mean[:, :, 0])
        & np.isfinite(sam_mean[:, :, 1])
    )
    if not np.any(valid):
        return math.nan, math.nan, 0
    w_error = candidate_w[valid] - sam_mean[:, :, 0][valid]
    rt_error = candidate_rt[valid] - sam_mean[:, :, 1][valid]
    return (
        float(np.sqrt(np.mean(w_error**2))),
        float(1000.0 * np.sqrt(np.mean(rt_error**2))),
        int(np.count_nonzero(valid)),
    )


def occupancy_comparison_figure():
    """Compare speed and equal-occupancy centers over all cached ARM/BOMEX data."""
    manifest = _manifest()
    case_names = [item["name"] for item in manifest["cases"]]
    labels = [name.upper() for name in case_names]
    figure = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Vertical-velocity center", "Total-water center"),
    )
    for scheme, label, color in (
        ("speed", "Crossing speed w", COLORS["up"]),
        ("uniform", "Occupancy (equal crossing)", COLORS["candidate"]),
    ):
        metrics = [_whole_case_center_metrics(name, scheme) for name in case_names]
        custom = np.array([[item[2]] for item in metrics], dtype=int)
        figure.add_trace(
            go.Bar(
                x=labels,
                y=[item[0] for item in metrics],
                name=label,
                legendgroup=scheme,
                marker_color=color,
                customdata=custom,
                hovertemplate="%{x}<br>RMSE=%{y:.3g} m/s<br>%{customdata[0]:,} scored planes<extra>%{fullData.name}</extra>",
            ),
            row=1,
            col=1,
        )
        figure.add_trace(
            go.Bar(
                x=labels,
                y=[item[1] for item in metrics],
                name=label,
                legendgroup=scheme,
                marker_color=color,
                customdata=custom,
                showlegend=False,
                hovertemplate="%{x}<br>RMSE=%{y:.3g} g/kg<br>%{customdata[0]:,} scored planes<extra>%{fullData.name}</extra>",
            ),
            row=1,
            col=2,
        )
    figure.update_yaxes(title_text="RMSE [m/s]", row=1, col=1)
    figure.update_yaxes(title_text="RMSE [g/kg]", row=1, col=2)
    figure.update_layout(
        template="plotly_dark",
        barmode="group",
        height=500,
        title="Whole-cache center test against the raw-SAM cloud-water-weighted population",
        legend={"orientation": "h", "y": -0.18},
        uirevision="pdf9-occupancy-comparison",
    )
    return figure


def _moment_preserving_mixture(snapshot, candidate, requested_fraction):
    packet_mean = candidate["mean"]
    packet_cov = candidate["covariance"]
    global_mean = snapshot.mean
    global_cov = snapshot.covariance
    fraction = float(np.clip(requested_fraction, 1.0e-4, 0.45))
    while fraction >= 1.0e-4:
        background_mean = (global_mean - fraction * packet_mean) / (1.0 - fraction)
        packet_delta = packet_mean - global_mean
        background_delta = background_mean - global_mean
        background_cov = (
            global_cov
            - fraction * (packet_cov + np.outer(packet_delta, packet_delta))
            - (1.0 - fraction) * np.outer(background_delta, background_delta)
        ) / (1.0 - fraction)
        if np.min(np.linalg.eigvalsh(0.5 * (background_cov + background_cov.T))) >= -1.0e-10:
            values, vectors = np.linalg.eigh(0.5 * (background_cov + background_cov.T))
            background_cov = vectors @ np.diag(np.maximum(values, 1.0e-14)) @ vectors.T
            return fraction, np.array((fraction, 1.0 - fraction)), np.vstack((packet_mean, background_mean)), np.stack((packet_cov, background_cov))
        fraction *= 0.72
    return 0.0, None, None, None


def closure_profile_figure(cube, meta, ti, fraction):
    """Compare raw-SAM and reconstructed cloud diagnostics through one column."""

    try:
        _step, raw_height, pressure = load_snapshot_pressure_profile(
            str(meta["sam_run_dir"]), int(round(float(cube["sam_time_seconds"][ti])))
        )
    except Exception as error:
        return _empty("Cloud-diagnostic reconstruction through the column", str(error))

    height = np.asarray(cube["height_m"], float)
    target = np.asarray(cube["sam_targets"][ti], float)
    predicted = np.full_like(target, np.nan, dtype=float)
    retained_fraction = np.full(len(height), np.nan)
    for di, selected_height in enumerate(height):
        candidate = _ensemble(cube, ti, di, OFFLINE_REFERENCE_WEIGHTING)
        if not candidate["available"]:
            continue
        mean = np.asarray(cube["sam_mean"][ti, di], float)
        covariance = np.asarray(cube["sam_covariance"][ti, di], float)
        raw_index = int(np.argmin(np.abs(raw_height - float(cube["sam_height_m"][di]))))
        chi = derive_chi_moments(
            mean[1],
            mean[2],
            float(pressure[raw_index]),
            covariance[1, 1],
            covariance[2, 2],
            covariance[1, 2],
            covariance[0, 1],
            covariance[0, 2],
            covariance[0, 0],
        )
        snapshot = SimpleNamespace(
            mean=mean,
            covariance=covariance,
            chi_mean=float(chi["mean_chi"]),
            chi_coeff_rt=float(chi["coef_rt"]),
            chi_coeff_thl=float(chi["coef_thl"]),
        )
        used, weights, means, covariances = _moment_preserving_mixture(
            snapshot, candidate, fraction
        )
        if used <= 0.0:
            continue
        predicted[di] = diagnose_cloud(weights, means, covariances, snapshot)[
            "diagnostics"
        ]
        retained_fraction[di] = used

    labels = ("w′rᶜ′", "w′²rᶜ′", "rₜ′rᶜ′", "θₗ′rᶜ′")
    figure = make_subplots(
        rows=1,
        cols=4,
        shared_yaxes=True,
        subplot_titles=labels,
    )
    for column, label in enumerate(labels, start=1):
        index = column - 1
        figure.add_trace(
            go.Scatter(
                x=target[:, index],
                y=height,
                mode="lines",
                line={"color": COLORS["sam"], "width": 3},
                name="Raw SAM target",
                legendgroup="sam",
                showlegend=column == 1,
                connectgaps=False,
            ),
            row=1,
            col=column,
        )
        figure.add_trace(
            go.Scatter(
                x=predicted[:, index],
                y=height,
                mode="lines+markers",
                line={"color": COLORS["candidate"], "width": 3},
                marker={"size": 4},
                name="Packet + residual Gaussian",
                legendgroup="reconstruction",
                showlegend=column == 1,
                customdata=retained_fraction,
                hovertemplate=(
                    "z=%{y:.0f} m<br>value=%{x:.4g}<br>retained ξ=%{customdata:.3f}"
                    "<extra>reconstruction</extra>"
                ),
                connectgaps=False,
            ),
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text=label, row=1, col=column)
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(
        template="plotly_dark",
        height=580,
        title=(
            "Cloud-diagnostic reconstruction through the column · "
            f"{cube['time_seconds'][ti] / 60.0:.0f} min"
        ),
        legend={"orientation": "h", "y": -0.18},
        uirevision="pdf9-closure-profile",
    )
    return figure


def continuity_figure(case_name, cube, scheme):
    w, rt, neff = _continuity_fields(case_name, scheme)
    rt = rt * 1000.0
    figure = make_subplots(rows=1, cols=3, subplot_titles=("Candidate rₜ", "Candidate w", "Effective parcel count"), shared_yaxes=True)
    x = cube["time_seconds"] / 60.0
    for col, field, title, colorscale in ((1, rt.T, "g/kg", "Viridis"), (2, w.T, "m/s", "RdBu"), (3, neff.T, "N_eff", "Purples")):
        figure.add_trace(go.Heatmap(x=x, y=cube["height_m"], z=field, colorscale=colorscale, colorbar={"title": title, "x": 0.29 + 0.35 * (col - 1), "len": 0.75}), row=1, col=col)
    figure.update_xaxes(title_text="LES time [min]")
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_layout(template="plotly_dark", height=560, title="Does the diagnosed transport population evolve continuously?", uirevision="pdf9-continuity")
    return figure


def recursive_candidate_figure(cube):
    """Show candidate provenance and when it feeds directional launches."""
    if "candidate_valid_combined" not in cube:
        return _empty("Next recursive G3 candidate", "Regenerate the PDF-9 cache")
    time = cube["time_seconds"] / 60.0
    height = cube["height_m"]
    valid = np.asarray(cube["candidate_valid_combined"], float) > 0.5
    upward_launch = np.asarray(cube["candidate_launch_from_g3_up"], float)
    downward_launch = np.asarray(cube["candidate_launch_from_g3_down"], float)
    count = np.asarray(cube["candidate_crossing_count_combined"], float)
    support = np.asarray(cube["candidate_weighted_support_combined"], float)
    distance = np.where(
        valid,
        np.asarray(cube["candidate_donor_distance_combined"], float),
        np.nan,
    )
    rt = np.where(
        valid, np.asarray(cube["candidate_rt_combined"], float) * 1000.0, np.nan
    )
    w = np.where(valid, np.asarray(cube["candidate_w_combined"], float), np.nan)
    probability_up = np.asarray(cube["candidate_branch_prob_up"], float)
    probability_down = np.asarray(cube["candidate_branch_prob_down"], float)
    figure = make_subplots(
        rows=2,
        cols=5,
        subplot_titles=(
            "Candidate rₜ",
            "Candidate w",
            "Raw crossings",
            "Branch/distance-weighted support",
            "Mean donor distance",
            "Prior G3 used for upward launch",
            "Prior G3 used for downward launch",
            "Positive-w branch probability",
            "Negative-w branch probability",
        ),
        shared_yaxes=True,
    )
    fields = (
        (1, 1, rt.T, "Viridis", "g/kg"),
        (1, 2, w.T, "Turbo", "m/s"),
        (1, 3, count.T, "Purples", "count"),
        (1, 4, support.T, "Purples", "weighted count"),
        (1, 5, distance.T, "Cividis", "m"),
        (
            2,
            1,
            upward_launch.T,
            ((0.0, "#111827"), (1.0, COLORS["recursive"])),
            "0/1",
        ),
        (
            2,
            2,
            downward_launch.T,
            ((0.0, "#111827"), (1.0, COLORS["down"])),
            "0/1",
        ),
        (2, 3, probability_up.T, "Blues", "P(w>0)"),
        (2, 4, probability_down.T, "Reds", "P(w<0)"),
    )
    for row, column, field, colorscale, title in fields:
        figure.add_trace(
            go.Heatmap(
                x=time,
                y=height,
                z=field,
                colorscale=colorscale,
                showscale=False,
                hovertemplate="time=%{x:.0f} min<br>z=%{y:.0f} m<br>value=%{z:.4g}<extra></extra>",
            ),
            row=row,
            col=column,
        )
    figure.update_xaxes(title_text="LES time [min]")
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_yaxes(title_text="Altitude [m]", row=2, col=1)
    figure.update_layout(
        template="plotly_dark",
        height=840,
        title=(
            "Distance-weighted bidirectional G3 and its conditional launch "
            "branches"
        ),
        uirevision="pdf9-recursive-candidate",
    )
    return figure


def downward_arrival_cap_figure(cube):
    """Show exactly where destination-local RMS clipping changes G3 pooling."""
    required = (
        "candidate_w_down_uncapped",
        "candidate_w_down",
        "candidate_down_cap_fraction",
        "candidate_destination_sigma_w",
    )
    if not all(name in cube for name in required):
        return _empty(
            "Destination-local downward-arrival cap",
            "Regenerate the PDF-9 cache to load the cap diagnostics",
        )
    time = np.asarray(cube["time_seconds"], float) / 60.0
    height = np.asarray(cube["height_m"], float)
    valid = np.asarray(cube["candidate_valid_down"], float) > 0.5
    raw = np.where(
        valid, np.asarray(cube["candidate_w_down_uncapped"], float), np.nan
    )
    capped = np.where(valid, np.asarray(cube["candidate_w_down"], float), np.nan)
    fraction = np.where(
        valid, np.asarray(cube["candidate_down_cap_fraction"], float), np.nan
    )
    threshold = -np.asarray(cube["candidate_destination_sigma_w"], float)
    sam_mean, _sam_covariance, sam_available = _sam_reference(
        cube, slice(None), "cloud_down"
    )
    sam_w = np.where(sam_available, sam_mean[:, :, 0], np.nan)
    finite_w = np.concatenate(
        tuple(field[np.isfinite(field)] for field in (raw, capped, sam_w))
    )
    bound = max(float(np.nanpercentile(np.abs(finite_w), 98.0)), 0.25) \
        if finite_w.size else 1.0
    figure = make_subplots(
        rows=2,
        cols=3,
        subplot_titles=(
            "Raw parcel-pool w",
            "Destination-capped w",
            "Raw SAM sinking-cloud w",
            "Applied center correction",
            "Support that hit the cap",
            "Destination limit −√w′²",
        ),
        shared_yaxes=True,
    )
    fields = (
        (1, 1, raw, "RdBu", -bound, bound, "m/s"),
        (1, 2, capped, "RdBu", -bound, bound, "m/s"),
        (1, 3, sam_w, "RdBu", -bound, bound, "m/s"),
        (2, 1, capped - raw, "Viridis", 0.0, None, "m/s"),
        (2, 2, fraction, "Magma", 0.0, 1.0, "fraction"),
        (2, 3, threshold, "RdBu", -bound, bound, "m/s"),
    )
    for row, column, field, colorscale, zmin, zmax, units in fields:
        figure.add_trace(
            go.Heatmap(
                x=time,
                y=height,
                z=np.asarray(field, float).T,
                colorscale=colorscale,
                zmin=zmin,
                zmax=zmax,
                showscale=False,
                hovertemplate=(
                    "time=%{x:.0f} min<br>z=%{y:.0f} m"
                    f"<br>value=%{{z:.4g}} {units}<extra></extra>"
                ),
            ),
            row=row,
            col=column,
        )
    figure.update_xaxes(title_text="LES time [min]")
    figure.update_yaxes(title_text="Altitude [m]", row=1, col=1)
    figure.update_yaxes(title_text="Altitude [m]", row=2, col=1)
    figure.update_layout(
        template="plotly_dark",
        height=820,
        title=(
            "Downward trajectories stay untouched; only their authority over "
            "the next G3 velocity moments is limited"
        ),
        uirevision="pdf9-downward-arrival-cap",
    )
    return figure


def _slider_marks(values, formatter):
    indices = np.unique(np.rint(np.linspace(0, len(values) - 1, min(7, len(values)))).astype(int))
    return {int(index): formatter(float(values[index])) for index in indices}


def _rail(manifest=None, cache_error: str = ""):
    cases = manifest.get("cases", []) if manifest else []
    first_name = cases[0]["name"] if cases else None
    if first_name is not None:
        _, cube = _case(first_name)
        times = cube["time_seconds"]
        heights = cube["height_m"]
    else:
        times = np.array((0.0,))
        heights = np.array((0.0,))
    return html.Aside(
        [
            html.Div([
                html.H3("Snapshot"),
                html.Label("Case"),
                dcc.Dropdown(id=CASE_ID, options=[{"label": item["name"].upper(), "value": item["name"]} for item in cases], value=first_name, clearable=False),
                html.Div([html.Span("LES time"), html.Strong(id=TIME_LABEL_ID)], className="coherent-control-heading"),
                dcc.Slider(id=TIME_ID, min=0, max=len(times)-1, step=1, value=len(times)//2, marks=_slider_marks(times, lambda value: f"{value/60:.0f}m"), updatemode="drag"),
                html.Div([html.Span("Altitude"), html.Strong(id=HEIGHT_LABEL_ID)], className="coherent-control-heading"),
                dcc.Slider(id=HEIGHT_ID, min=0, max=len(heights)-1, step=1, value=min(len(heights)//3, len(heights)-1), marks=_slider_marks(heights, lambda value: f"{value:.0f}m"), updatemode="drag"),
            ], className="coherent-rail-section coherent-rail-controls"),
            html.Div([html.H3("Selected destination"), html.Div(id=FACTS_ID)], className="coherent-rail-section"),
        ], className="coherent-control-rail"
    )


def _generation_controls(cache_error: str = ""):
    """Keep the expensive regeneration action visible without persistent log text."""

    return html.Section(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.H3("PDF-9 diagnostics data"),
                            html.P(
                                "Run fresh ARM and BOMEX PDF-9 cases, then rebuild the Misc cache.",
                                className="coherent-control-note",
                            ),
                        ]
                    ),
                    html.Button(
                        "Generate / refresh PDF-9 data",
                        id=REFRESH_ID,
                        n_clicks=0,
                        className="coherent-action-button",
                    ),
                ],
                className="pdf9-generation-heading",
            ),
            html.P(
                cache_error,
                id=CACHE_STATUS_ID,
                className="coherent-control-note pdf9-generation-status",
                style={} if cache_error else {"display": "none"},
            ),
        ],
        className="pdf9-generation-bar",
    )


def build_layout():
    try:
        manifest = _manifest()
        _validate_cache_files(manifest)
    except Exception as error:
        manifest = None
        cache_error = str(error)
    else:
        cache_error = ""
    sections = html.Div([
        report_section("1. Candidate Gaussians over raw SAM", html.Div([html.Div([html.Label("Horizontal variable"), dcc.Dropdown(id=OVERLAY_VARIABLE_ID, options=[{"label": value["label"], "value": key} for key, value in OVERLAY_VARIABLES.items()], value="rt", clearable=False)], className="pdf9-overlay-selector"), dcc.Graph(id=OVERLAY_ID, config={"displaylogo": False})]), intro="Choose w–rₜ, w–χ, or w–r_c. The rₜ and χ views show component contours; r_c is nonlinear, so its view shows the raw SAM structure and PDF component centers.", collapsible=True, open_by_default=True),
        report_section("2. PDF component weights through the column", dcc.Graph(id=COMPONENT_WEIGHTS_ID, config={"displaylogo": False}), intro="The three lines are the actual PDF mixture weights: G1 and G2 share the residual mass, while G3 is the recursive transport component.", collapsible=True, open_by_default=True),
        report_section("3. Follow every parcel track", dcc.Graph(id=TRACKS_ID, config={"displaylogo": False}), intro="Faint cyan and rose curves are upward and downward launch parcels. Where pooled G3 exists, the rising and sinking solvers use its positive- and negative-w conditional means capped at one local vertical-velocity RMS; rₜ and θₗ are conditioned at that capped velocity. Lime follows the pooled G3 through the column; gold and dotted cyan are evaluation-only SAM references.", collapsible=True, open_by_default=True),
        report_section("4. Compare mixing reach and source support", dcc.Graph(id=REACH_ID, config={"displaylogo": False}), intro="The left plot overlays upward and downward reach from every launch. The other two show which source heights support each destination and the corresponding effective parcel count.", collapsible=True, open_by_default=True),
        report_section("5. Reconstruct cloud diagnostics through the column", dcc.Graph(id=CLOSURE_ID, config={"displaylogo": False}), intro="At the selected time, each panel overlays a raw-SAM cloud diagnostic and the packet-plus-residual reconstruction across altitude. This is a capacity test, not a diagnosed PDF weight.", collapsible=True, open_by_default=True),
        report_section("6. Audit the next recursive G3 candidate", dcc.Graph(id=RECURSIVE_ID, config={"displaylogo": False}), intro="Every incoming crossing is weighted by its directional Gaussian probability mass times exp(-mu times donor distance) in the next G3 mean and full covariance. Raw crossing count, weighted support, and mean donor distance are shown separately. On the next PDF call, w′rₜ′ determines how much probability this geometry may receive; any unrealizable residual covariance backs that weight down.", collapsible=True, open_by_default=True),
        report_section("7. Audit the destination-velocity cap", dcc.Graph(id=CAP_ID, config={"displaylogo": False}), intro="The trajectory ledger and mixing reach are unchanged. For pooling only, each sinking arrival contributes max(wparcel, −sqrt(w′²destination)) to the next G3 velocity mean and covariance. The raw value, capped value, cap frequency, and held-out SAM sinking-cloud center remain visible here.", collapsible=True, open_by_default=True),
    ], className="coherent-results-pane")
    return html.Article([
        dcc.Store(id=CACHE_REVISION_ID, data=0, storage_type="memory"),
        dcc.Interval(
            id=GENERATION_INTERVAL_ID,
            interval=1_000,
            n_intervals=0,
            disabled=True,
        ),
        report_header("PDF-9 parcel-ensemble laboratory", "Turn the new launch-by-destination parcel ledger into a candidate transport Gaussian, then test its placement, shape, support, continuity, and cloud-moment capacity against raw ARM and BOMEX SAM.", eyebrow="Developing closure method", badges=("PDF type 9", "LES-advance", "SAM scoring only")),
        _generation_controls(cache_error),
        html.Div([sections, html.Div(className="coherent-pane-divider"), _rail(manifest, cache_error)], className="coherent-workspace"),
    ], className="notes-report coherent-report pdf9-report")


def register_callbacks(app):
    @app.callback(
        Output(CASE_ID, "options"),
        Output(CASE_ID, "value"),
        Output(CACHE_STATUS_ID, "children"),
        Output(CACHE_STATUS_ID, "style"),
        Output(CACHE_REVISION_ID, "data"),
        Output(GENERATION_INTERVAL_ID, "disabled"),
        Input(REFRESH_ID, "n_clicks"),
        Input(GENERATION_INTERVAL_ID, "n_intervals"),
        State(CASE_ID, "value"),
        State(CACHE_REVISION_ID, "data"),
        prevent_initial_call=True,
    )
    def generate_or_poll(_clicks, _ticks, selected_case, revision):
        if ctx.triggered_id == REFRESH_ID:
            return (
                no_update,
                no_update,
                start_pdf9_generation(),
                {},
                no_update,
                False,
            )

        running, succeeded, message = job_status("pdf9_parcel_ensemble")
        if running:
            return no_update, no_update, message, {}, no_update, False
        if not succeeded:
            return no_update, no_update, message or no_update, {}, no_update, True
        _clear_data_caches()
        try:
            manifest = _manifest()
        except Exception as error:
            return (
                [],
                None,
                f"{message}\n\nCache validation failed: {error}",
                {},
                int(revision or 0) + 1,
                True,
            )
        cases = [item["name"] for item in manifest["cases"]]
        selected = selected_case if selected_case in cases else cases[0]
        return (
            [{"label": name.upper(), "value": name} for name in cases],
            selected,
            "",
            {"display": "none"},
            int(revision or 0) + 1,
            True,
        )

    @app.callback(
        Output(TIME_ID, "max"), Output(TIME_ID, "marks"), Output(TIME_ID, "value"),
        Output(HEIGHT_ID, "max"), Output(HEIGHT_ID, "marks"), Output(HEIGHT_ID, "value"),
        Input(CASE_ID, "value"), Input(CACHE_REVISION_ID, "data"),
    )
    def change_case(case_name, _revision):
        if not case_name:
            return 0, {}, 0, 0, {}, 0
        _, cube = _case(case_name)
        return (
            len(cube["time_seconds"])-1, _slider_marks(cube["time_seconds"], lambda value: f"{value/60:.0f}m"), len(cube["time_seconds"])//2,
            len(cube["height_m"])-1, _slider_marks(cube["height_m"], lambda value: f"{value:.0f}m"), min(len(cube["height_m"])//3, len(cube["height_m"])-1),
        )

    @app.callback(
        Output(TIME_LABEL_ID, "children"), Output(HEIGHT_LABEL_ID, "children"), Output(FACTS_ID, "children"),
        Output(OVERLAY_ID, "figure"), Output(COMPONENT_WEIGHTS_ID, "figure"), Output(TRACKS_ID, "figure"), Output(REACH_ID, "figure"), Output(CLOSURE_ID, "figure"), Output(RECURSIVE_ID, "figure"), Output(CAP_ID, "figure"),
        Input(CASE_ID, "value"), Input(TIME_ID, "value"), Input(HEIGHT_ID, "value"), Input(OVERLAY_VARIABLE_ID, "value"), Input(CACHE_REVISION_ID, "data"),
    )
    def render(case_name, time_index, height_index, overlay_variable, _revision):
        scheme = OFFLINE_REFERENCE_WEIGHTING
        fraction = CLOSURE_PACKET_FRACTION
        if not case_name:
            figure = _empty("PDF-9 parcel laboratory", "Load a PDF-9 cache")
            return (
                "—", "—",
                empty_state("No cache loaded", "Generate the cache and click Reload cache"),
                *(figure for _ in range(7)),
            )
        try:
            meta, cube = _case(case_name)
        except Exception as error:
            figure = _empty("PDF-9 parcel laboratory", str(error))
            return (
                "—", "—",
                empty_state("Cache unavailable", str(error)),
                *(figure for _ in range(7)),
            )
        ti = int(np.clip(time_index or 0, 0, len(cube["time_seconds"])-1))
        di = int(np.clip(height_index or 0, 0, len(cube["height_m"])-1))
        overlay, _snapshot, candidate = overlay_figure(
            cube, meta, ti, di, scheme, overlay_variable or "rt"
        )
        downward_reduction = _downward_ensemble(cube, ti, di)
        downward = _stored_candidate(cube, ti, di, "down")
        sam_mean, _, available = _sam_reference(cube, ti)
        if candidate["available"] and available[di]:
            error_w = candidate["mean"][0] - sam_mean[di, 0]
            error_rt = (candidate["mean"][1] - sam_mean[di, 1]) * 1000.0
        else:
            error_w = error_rt = math.nan
        crossing_count = int(np.count_nonzero(cube["parcel_status"][ti, di] > 0))
        down_count = int(np.count_nonzero(cube["down_parcel_status"][ti, di] > 0))
        recursive = _recursive_candidate(cube, ti, di)
        active_g3 = _active_g3(cube, ti, di)
        recursive_count = _finite_int(
            cube["candidate_crossing_count_combined"][ti, di]
        )
        recursive_support = _finite_float(
            cube["candidate_weighted_support_combined"][ti, di]
        )
        recursive_distance = _finite_float(
            cube["candidate_donor_distance_combined"][ti, di], math.nan
        )
        downward_support = _finite_float(
            cube["candidate_weighted_support_down"][ti, di]
        )
        downward_candidate_count = _finite_int(
            cube["candidate_crossing_count_down"][ti, di]
        )
        recursive_launch = bool(
            cube["candidate_launch_from_g3_up"][ti, di] > 0.5
        )
        recursive_down_launch = bool(
            cube["candidate_launch_from_g3_down"][ti, di] > 0.5
        )
        down_mean, _, down_available = _sam_reference(cube, ti, "cloud_down")
        if downward["available"] and down_available[di]:
            down_error_w = downward["mean"][0] - down_mean[di, 0]
            down_error_rt = (downward["mean"][1] - down_mean[di, 1]) * 1000.0
        else:
            down_error_w = down_error_rt = math.nan
        raw_down_w = _finite_float(
            cube["candidate_w_down_uncapped"][ti, di], math.nan
        )
        capped_down_w = _finite_float(cube["candidate_w_down"][ti, di], math.nan)
        cap_fraction = _finite_float(cube["candidate_down_cap_fraction"][ti, di])
        destination_sigma = _finite_float(
            cube["candidate_destination_sigma_w"][ti, di], math.nan
        )
        facts = fact_grid([
            {"label":"Upward crossings", "value":str(crossing_count), "detail":f"N_eff {candidate['neff']:.2f}; {candidate['count']} weighted"},
            {"label":"Downward crossings", "value":str(down_count), "detail":f"Fortran support {downward_support:.2f} from {downward_candidate_count} arrivals", "tone":"good" if downward_candidate_count >= 3 else "warning"},
            {"label":"Upward center error", "value":f"{error_w:+.2f} m/s", "detail":f"rₜ {error_rt:+.2f} g/kg"},
            {"label":"Downward center error", "value":f"{down_error_w:+.2f} m/s", "detail":f"rₜ {down_error_rt:+.2f} g/kg vs sinking cloud"},
            {"label":"Raw SAM cloud", "value":f"{cube['sam_cloud_fraction'][ti,di]*100:.2f}%", "detail":"evaluation only"},
            {"label":"Up/down reach", "value":f"{cube['lscale_up'][ti,di]:.0f} / {cube['lscale_down'][ti,di]:.0f} m", "detail":"grid-mean parcel launch"},
            {"label":"Up/down crossing sum", "value":f"{np.sum(candidate['weight']):.3g} / {np.sum(downward_reduction['weight']):.3g}", "detail":"reduction weight—not ξ"},
            {"label":"Next G3 candidate", "value":"valid" if recursive["available"] else "absent", "detail":f"support {recursive_support:.2f} from {recursive_count} crossings; donor distance {_distance_label(recursive_distance)}", "tone":"good" if recursive["available"] else "warning"},
            {"label":"Active G3 weight", "value":f"{active_g3['weight']:.4f}", "detail":"analytic share of w′rₜ′; residual moments removed exactly", "tone":"good" if active_g3["available"] else "warning"},
            {"label":"G3 entrainment blend", "value":f"{_finite_float(cube['entrain_g3_weight'][ti,di]):.4f}", "detail":f"Δrₜ up/down {(_finite_float(cube['entrain_rt_up'][ti,di])-_finite_float(cube['mean_rt'][ti,di]))*1000.0:+.3f} / {(_finite_float(cube['entrain_rt_down'][ti,di])-_finite_float(cube['mean_rt'][ti,di]))*1000.0:+.3f} g/kg"},
            {"label":"Upward launch source", "value":"pooled G3" if recursive_launch else "current moments", "detail":"positive-w conditional branch"},
            {"label":"Downward launch source", "value":"pooled G3" if recursive_down_launch else "current moments", "detail":"negative-w conditional branch of the same G3"},
            {"label":"Downward arrival cap", "value":f"{raw_down_w:.2f} → {capped_down_w:.2f} m/s", "detail":f"{cap_fraction*100.0:.0f}% of support capped at −{destination_sigma:.2f} m/s", "tone":"good" if cap_fraction > 0.0 else "neutral"},
        ])
        return (
            f"{cube['time_seconds'][ti]/60:.0f} min", f"{cube['height_m'][di]:.0f} m", facts,
            overlay, component_weights_figure(cube, ti), tracks_figure(cube, ti, di), reach_and_support_figure(cube, ti, di, scheme), closure_profile_figure(cube, meta, ti, float(fraction)), recursive_candidate_figure(cube), downward_arrival_cap_figure(cube),
        )


REPORT = register_report(ReportSpec(
    slug="pdf9-parcel-ensemble",
    title="PDF-9 parcel-ensemble laboratory",
    summary="Use the launch-by-destination parcel ledger to diagnose a transport Gaussian and audit its support, geometry, continuity, and cloud-moment capacity.",
    category="PDF development",
    updated="2026-07-19",
    tags=("PDF 9", "parcel tracking", "SAM 3-D", "Gaussian 3"),
    order=37,
    build_layout=build_layout,
    register_callbacks=register_callbacks,
))
