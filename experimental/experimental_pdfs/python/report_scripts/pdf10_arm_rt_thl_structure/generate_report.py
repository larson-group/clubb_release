#!/usr/bin/env python3
"""Render the static ARM r_t--theta_l PDF-10 investigation report.

The figures are deliberately rendered here and copied into the report bundle;
the Dash Reports tab only serves the finished artifacts.  This script uses the
current full-LES-advance PDF-10 ARM output plus the matching raw 3-D SAM
planes.  Re-run only when intentionally revising the investigation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import xarray as xr


REPORT_ID = "pdf10-arm-rt-thl-structure"
REPORT_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPORT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_pdf10_rt_thl_report_matplotlib")

from dash_app.reports_tab.publisher import ReportBuilder
from utilities.sam_3d_reference import DEFAULT_SAM_RUN, load_snapshot


STATS_PATH = REPO_ROOT / "output" / "arm_stats.nc"
LES_PATH = REPO_ROOT / "output" / "arm_les_advance.nc"
ARM_INITIAL_SECONDS = 41_400.0
SELECTED_HOURS = (15.5, 18.0, 20.5)
SELECTED_HEIGHTS = (820.0, 1_500.0, 2_100.0)
MASS_LEVELS = (0.50, 0.80, 0.995)
COMPONENT_COLORS = ("#ff5c67", "#28d7f5")


def _git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip() + " (working tree)"
    except (OSError, subprocess.CalledProcessError):
        return "working tree"


def _scalar(data_array, time_index: int, height_index: int) -> float:
    values = data_array.isel(time=time_index, **{data_array.dims[1]: height_index})
    for dim in tuple(values.dims):
        values = values.isel({dim: 0})
    return float(values.values)


def _nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def _mass_levels(probability: np.ndarray) -> list[float]:
    flat = np.asarray(probability, dtype=float).ravel()
    flat = flat[np.isfinite(flat) & (flat > 0.0)]
    if flat.size == 0:
        return []
    ranked = np.sort(flat)[::-1]
    cumulative = np.cumsum(ranked)
    total = cumulative[-1]
    return [float(ranked[min(np.searchsorted(cumulative, mass * total), ranked.size - 1)]) for mass in MASS_LEVELS]


def _histogram(first: np.ndarray, second: np.ndarray, bins: int = 74):
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    x_low, x_high = np.percentile(first, (0.15, 99.85))
    y_low, y_high = np.percentile(second, (0.15, 99.85))
    x_pad = max(0.055 * (x_high - x_low), 1.0e-5)
    y_pad = max(0.055 * (y_high - y_low), 1.0e-5)
    x_edges = np.linspace(x_low - x_pad, x_high + x_pad, bins + 1)
    y_edges = np.linspace(y_low - y_pad, y_high + y_pad, bins + 1)
    probability = np.histogram2d(first, second, bins=(x_edges, y_edges))[0].T / first.size
    return x_edges, y_edges, probability


def _ellipse(center: tuple[float, float], covariance: np.ndarray, radius: float = 2.0):
    covariance = 0.5 * (np.asarray(covariance, dtype=float) + np.asarray(covariance, dtype=float).T)
    values, vectors = np.linalg.eigh(covariance)
    values = np.maximum(values, 0.0)
    theta = np.linspace(0.0, 2.0 * np.pi, 180)
    points = vectors @ (radius * np.sqrt(values)[:, None] * np.stack((np.cos(theta), np.sin(theta))))
    return center[0] + points[0], center[1] + points[1]


def _component_projection(stats, time_index: int, height_index: int, first: str, second: str, component: int):
    suffix = str(component)
    means = {
        "w": _scalar(stats[f"w_{suffix}"], time_index, height_index),
        "rt": _scalar(stats[f"rt_{suffix}"], time_index, height_index) * 1_000.0,
        "thl": _scalar(stats[f"thl_{suffix}"], time_index, height_index),
    }
    variances = {
        "w": _scalar(stats[f"varnce_w_{suffix}"], time_index, height_index),
        "rt": _scalar(stats[f"varnce_rt_{suffix}"], time_index, height_index) * 1.0e6,
        "thl": _scalar(stats[f"varnce_thl_{suffix}"], time_index, height_index),
    }
    corr_name = {
        frozenset(("w", "rt")): f"corr_w_rt_{suffix}",
        frozenset(("w", "thl")): f"corr_w_thl_{suffix}",
        frozenset(("rt", "thl")): f"corr_rt_thl_{suffix}",
    }[frozenset((first, second))]
    correlation = _scalar(stats[corr_name], time_index, height_index)
    covariance = correlation * np.sqrt(max(variances[first], 0.0) * max(variances[second], 0.0))
    return (means[first], means[second]), np.array(
        ((variances[first], covariance), (covariance, variances[second])), dtype=float
    ), correlation


def _raw_axes(snapshot, first: str, second: str):
    data = {
        "w": np.asarray(snapshot.samples[:, 0], dtype=float),
        "rt": np.asarray(snapshot.samples[:, 1], dtype=float) * 1_000.0,
        "thl": np.asarray(snapshot.samples[:, 2], dtype=float),
    }
    return data[first], data[second]


def _draw_plane(axis, first, second, snapshot, stats, time_index, height_index, *, title, show_labels=True):
    from matplotlib.colors import LogNorm

    x, y = _raw_axes(snapshot, first, second)
    x_edges, y_edges, probability = _histogram(x, y)
    nonzero = probability[probability > 0.0]
    axis.set_facecolor("#111827")
    axis.pcolormesh(
        x_edges,
        y_edges,
        probability,
        shading="auto",
        cmap="magma",
        norm=LogNorm(vmin=max(float(np.percentile(nonzero, 5)), 1.0e-8), vmax=float(np.max(nonzero))),
    )
    levels = _mass_levels(probability)
    if levels:
        axis.contour(
            0.5 * (x_edges[:-1] + x_edges[1:]),
            0.5 * (y_edges[:-1] + y_edges[1:]),
            probability,
            levels=sorted(levels),
            colors="#facc15",
            linewidths=(1.0, 1.4, 1.9),
        )
    for component, color in enumerate(COMPONENT_COLORS, start=1):
        mean, covariance, correlation = _component_projection(
            stats, time_index, height_index, first, second, component
        )
        if np.isfinite(covariance).all() and np.isfinite(mean).all():
            ellipse_x, ellipse_y = _ellipse(mean, covariance)
            axis.plot(ellipse_x, ellipse_y, color=color, linewidth=2.0, linestyle="--")
            axis.scatter(*mean, color=color, marker="x", s=38, linewidth=1.7, zorder=5)
    axis.set_title(title, color="#edf2f7", fontsize=9.5, pad=7)
    axis.tick_params(colors="#cbd5e1", labelsize=8)
    for spine in axis.spines.values():
        spine.set_color("#40516a")
    if show_labels:
        labels = {"w": "w [m/s]", "rt": "rₜ [g/kg]", "thl": "θₗ [K]"}
        axis.set_xlabel(labels[first], color="#e5e7eb", fontsize=8.5)
        axis.set_ylabel(labels[second], color="#e5e7eb", fontsize=8.5)
    return {
        "raw_correlation": float(np.corrcoef(x, y)[0, 1]),
        "component_correlations": [
            _component_projection(stats, time_index, height_index, first, second, component)[2]
            for component in (1, 2)
        ],
    }


def _model_time(hours: float) -> float:
    return float(hours) * 3600.0


def _raw_snapshot_for_time(model_seconds: float, height: float):
    return load_snapshot(DEFAULT_SAM_RUN, int(round(model_seconds - ARM_INITIAL_SECONDS)), height)


def _render_rt_thl_mosaic(work_dir: Path, stats) -> list[dict]:
    from matplotlib import pyplot as plt

    figure, axes = plt.subplots(3, 3, figsize=(15.5, 13.6), dpi=170, facecolor="#0d1525")
    rows = []
    for row, height in enumerate(SELECTED_HEIGHTS):
        for column, hours in enumerate(SELECTED_HOURS):
            model_seconds = _model_time(hours)
            time_index = _nearest_index(stats.time.values, model_seconds)
            height_index = _nearest_index(stats.zt.values, height)
            snapshot = _raw_snapshot_for_time(model_seconds, height)
            info = _draw_plane(
                axes[row, column], "rt", "thl", snapshot, stats, time_index, height_index,
                title=(f"{hours:.1f} h / {snapshot.elapsed_minutes:.0f} raw min\n"
                       f"z={snapshot.height_m:.0f} m, raw ρ={np.corrcoef(snapshot.samples[:,1], snapshot.samples[:,2])[0,1]:+.3f}"),
            )
            rows.append({
                "model_time_seconds": model_seconds,
                "model_time_hours": hours,
                "raw_elapsed_minutes": snapshot.elapsed_minutes,
                "height_m": snapshot.height_m,
                "sample_count": snapshot.sample_count,
                "raw_corr_rt_thl": info["raw_correlation"],
                "pdf10_component_corr_rt_thl": info["component_correlations"],
                "source_file": snapshot.source_file,
            })
    figure.suptitle(
        "ARM raw-SAM rₜ–θₗ planes with PDF-10 components (full LES advance)",
        color="#f8fafc", fontsize=16, y=0.995,
    )
    figure.text(
        0.5, 0.012,
        "Magma: raw-SAM bin probability. Gold: raw 50/80/99.5% highest-density contours. "
        "Dashed red/cyan: PDF-10 Gaussian 1/2 2σ shapes; crosses are their means.",
        ha="center", color="#b6c5d8", fontsize=9.5,
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.97))
    output = work_dir / "rt-thl-raw-sam-pdf10-mosaic.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    return rows


def _render_projection_comparison(work_dir: Path, stats) -> dict:
    from matplotlib import pyplot as plt

    hours = 20.5
    model_seconds = _model_time(hours)
    time_index = _nearest_index(stats.time.values, model_seconds)
    figure, axes = plt.subplots(3, 3, figsize=(15.5, 13.6), dpi=170, facecolor="#0d1525")
    planes = (("rt", "thl", "rₜ–θₗ"), ("w", "rt", "w–rₜ"), ("w", "thl", "w–θₗ"))
    result = {"model_time_hours": hours, "model_time_seconds": model_seconds, "planes": {}}
    for column, height in enumerate(SELECTED_HEIGHTS):
        height_index = _nearest_index(stats.zt.values, height)
        snapshot = _raw_snapshot_for_time(model_seconds, height)
        for row, (first, second, label) in enumerate(planes):
            info = _draw_plane(
                axes[row, column], first, second, snapshot, stats, time_index, height_index,
                title=f"{label}: z={snapshot.height_m:.0f} m, raw ρ={np.corrcoef(*_raw_axes(snapshot, first, second))[0,1]:+.3f}",
            )
            result["planes"].setdefault(label, []).append({
                "height_m": snapshot.height_m,
                "raw_correlation": info["raw_correlation"],
                "pdf10_component_correlations": info["component_correlations"],
            })
    figure.suptitle(
        "Why rₜ–θₗ looks different: matched raw-SAM projections at 20.5 h", color="#f8fafc", fontsize=16, y=0.995
    )
    figure.text(
        0.5, 0.012,
        "Same rendering convention in every panel.  The rₜ–θₗ row is persistently narrow and tilted; "
        "the two w planes carry substantially less coherent orientation.",
        ha="center", color="#b6c5d8", fontsize=9.5,
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.97))
    output = work_dir / "rt-thl-vs-w-projections.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    return result


def _render_forced_moment_context(work_dir: Path, les) -> dict:
    from matplotlib import pyplot as plt

    def field(name):
        values = np.asarray(les[name].values, dtype=float)
        return values[:, :, 0, 0]

    rt2, thl2, w2 = field("rtp2"), field("thlp2"), field("wp2")
    rtthl, wrt, wthl = field("rtpthlp"), field("wprtp"), field("wpthlp")
    def correlation(covariance, left, right):
        answer = covariance / np.sqrt(np.maximum(left * right, 1.0e-300))
        return np.where(np.isfinite(answer), np.clip(answer, -1.0, 1.0), np.nan)
    fields = (
        (correlation(rtthl, rt2, thl2), "ρ(rₜ, θₗ)"),
        (correlation(wrt, w2, rt2), "ρ(w, rₜ)"),
        (correlation(wthl, w2, thl2), "ρ(w, θₗ)"),
        (np.log10(np.maximum(rt2 * 1.0e6, 1.0e-14)), "log₁₀ var(rₜ) [(g/kg)²]"),
    )
    time_hours = np.asarray(les.time.values, dtype=float) / 3600.0
    z = np.asarray(les.z.values, dtype=float)
    figure, axes = plt.subplots(2, 2, figsize=(15.2, 10.6), dpi=170, facecolor="#0d1525", sharex=True, sharey=True)
    for axis, (values, title) in zip(axes.flat, fields):
        axis.set_facecolor("#111827")
        if title.startswith("ρ"):
            image = axis.pcolormesh(time_hours, z, values.T, shading="auto", cmap="coolwarm", vmin=-1.0, vmax=1.0)
        else:
            image = axis.pcolormesh(time_hours, z, values.T, shading="auto", cmap="viridis")
        figure.colorbar(image, ax=axis, pad=0.015, shrink=0.87).ax.tick_params(colors="#cbd5e1", labelsize=8)
        axis.set_title(title, color="#edf2f7", fontsize=11)
        axis.tick_params(colors="#cbd5e1", labelsize=8)
        for spine in axis.spines.values():
            spine.set_color("#40516a")
        for hours in SELECTED_HOURS:
            axis.axvline(hours, color="#facc15", alpha=0.62, linewidth=1.0, linestyle="--")
    for axis in axes[-1, :]:
        axis.set_xlabel("CLUBB/SAM absolute clock [h]", color="#e5e7eb")
    for axis in axes[:, 0]:
        axis.set_ylabel("height [m]", color="#e5e7eb")
    figure.suptitle(
        "ARM LES-forced covariance context supplied to PDF-10", color="#f8fafc", fontsize=16, y=0.99
    )
    figure.text(
        0.5, 0.012,
        "Dashed lines mark the static-plane examples.  Near-unit correlations above the well-variance-supported cloud layer must be read with the rₜ-variance panel; tiny variances make any correlation fragile.",
        ha="center", color="#b6c5d8", fontsize=9.5,
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.96))
    output = work_dir / "arm-les-forced-covariance-context.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    valid = np.isfinite(fields[0][0]) & (rt2 > 1.0e-14) & (thl2 > 1.0e-8)
    return {
        "corr_rt_thl_percentiles": [float(value) for value in np.nanpercentile(fields[0][0][valid], (1, 5, 25, 50, 75, 95, 99))],
        "valid_corr_rt_thl_count": int(np.count_nonzero(valid)),
        "total_states": int(valid.size),
    }


def _publish(work_dir: Path, selection: list[dict], projections: dict, context: dict):
    source = {
        "stats_file": str(STATS_PATH),
        "les_advance_file": str(LES_PATH),
        "raw_sam_run": str(DEFAULT_SAM_RUN),
        "arm_model_time_initial_seconds": ARM_INITIAL_SECONDS,
        "raw_time_mapping": "raw elapsed seconds = CLUBB absolute seconds - 41400",
        "selected_model_hours": list(SELECTED_HOURS),
        "selected_heights_m": list(SELECTED_HEIGHTS),
        "raw_probability_masses": list(MASS_LEVELS),
        "pdf_component_outline": "2 standard deviations of each bivariate component marginal",
    }
    data_path = work_dir / "rt-thl-structure-data.json"
    data_path.write_text(json.dumps({"source": source, "selected_planes": selection, "projection_comparison": projections, "context": context}, indent=2) + "\n", encoding="utf-8")
    report = ReportBuilder(
        REPORT_ID,
        "PDF-10 in ARM: what the rₜ–θₗ contours are telling us",
        summary="Matched raw-SAM planes and a full-LES-advance PDF-10 reconstruction show a persistent, nearly one-dimensional thermodynamic rₜ–θₗ direction that is much more coherent than the corresponding w relationships.",
        tags=("PDF-10", "ARM", "SAM", "LES advance", "r_t", "theta_l", "trivariate PDF"),
        replace=True,
        source_revision=_git_revision(),
    )
    report.heading("Answer at a glance")
    report.callout(
        "The useful signal is real — and it exposes a component-allocation weakness.",
        "Across the well-variance-supported ARM cloud layer, raw SAM and the LES-forced covariance both show a strongly negative r_t–theta_l relationship.  In several matched planes, PDF-10 Gaussian 2 follows that narrow tilt, while Gaussian 1 is much broader or displaced from the high-probability raw support.  Because all relevant grid moments are LES-forced here, that is a representation/allocation issue inside the closure rather than moment drift.  The direction is valuable for component shape and cloud diagnosis, but it is not evidence by itself that component centers should separate along that axis.",
        tone="warning",
    )
    report.metrics([
        ("Raw planes", "9 + 9", "nine r_t–theta_l examples plus nine matched w comparisons"),
        ("Forced inputs", "complete", "means, covariance matrix, wp3, rtp3, and thlp3 supplied from ARM LES"),
        ("r_t–theta_l median ρ", f"{context['corr_rt_thl_percentiles'][3]:+.2f}", "all variance-supported time-height states"),
        ("Interpretation", "shape axis", "retain as component covariance; do not automatically spend center separation on it"),
    ])
    report.heading("Matched rₜ–θₗ snapshots")
    report.paragraph(
        "Each panel uses an actual horizontal raw-SAM plane at the corresponding ARM time and height.  Gold contours enclose 50%, 80%, and 99.5% of the binned raw probability; PDF-10’s two component marginals are overlaid without fitting them to the plane.  Full LES advance makes this a representation test: the grid moments are supplied from SAM, while the component allocation and cloud diagnosis remain the closure’s work.  The repeated broad/displaced red ellipse is therefore a useful diagnostic: matching a total covariance does not guarantee that each component occupies plausible thermodynamic support."
    )
    report.figure(work_dir / "rt-thl-raw-sam-pdf10-mosaic.png", caption="Raw SAM r_t–theta_l probability geometry with the PDF-10 two-component representation.  The long negative tilt is persistent through the active cloud layer, while the raw distribution can still have shoulders, curvature, or unequal tails that two ellipses cannot reproduce exactly.")
    report.heading("Why this plane differs from the w projections")
    report.paragraph(
        "At 20.5 h, the same raw points are projected three ways.  The thermodynamic pair behaves like a compact mixing/saturation direction: moist anomalies and cool theta_l anomalies reinforce the saturation-excess coordinate.  The w projections retain plume information, but their orientation is more intermittent because vertical velocity has broader, reversible plume dynamics and must also carry the prescribed w variance."
    )
    report.figure(work_dir / "rt-thl-vs-w-projections.png", caption="Same raw SAM planes, time, heights, probability contours, and PDF-10 outlines.  The top r_t–theta_l row is visibly more coherent than the two w-based rows.")
    report.heading("Moment context and the important caveat")
    report.paragraph(
        "The input covariance supports the visual impression: corr(r_t,theta_l) is strongly negative through much of the cloud layer.  However, a correlation coefficient becomes uninformative when either variance is tiny.  The lower-right panel is included specifically to prevent interpreting noisy upper-level near-unit correlations as a robust thermodynamic constraint."
    )
    report.figure(work_dir / "arm-les-forced-covariance-context.png", caption="LES-advanced fields supplied to PDF-10.  The r_t–theta_l covariance is a widespread cloud-layer feature; compare it with the r_t variance before giving an extreme correlation physical weight.")
    report.heading("How PDF-10 can use this without becoming complicated")
    report.paragraph(
        "PDF-10 already preserves the supplied r_t–theta_l covariance through its component covariance matrices, so the immediate opportunity is not to add another predicted grid moment.  It is to use the known thermodynamic direction more deliberately when allocating residual covariance and width contrast between G1 and G2.  This matters strongly for cloud because the linearized saturation coordinate is chi' = c_rt r_t' - c_thl theta_l'.  A negative covariance therefore increases chi variance: the two terms reinforce rather than cancel."
    )
    report.equation(r"\sigma_\chi^2 = c_{rt}^2\sigma_{rt}^2 + c_{\theta_l}^2\sigma_{\theta_l}^2 - 2c_{rt}c_{\theta_l}\operatorname{cov}(r_t,\theta_l)", caption="Negative cov(r_t, theta_l) broadens the saturation-excess distribution when both coefficients are positive.")
    report.paragraph(
        "A conservative next closure experiment would first protect the observed narrow thermodynamic direction in both residual component covariances, while using the available skewness information to decide only modest center separation along a whitened thermodynamic coordinate.  That directly targets the broad/displaced-G1 symptom without confusing a narrow tilted ellipse with proof of two distinct thermodynamic populations.  A more ambitious center-separation rule should first be tested against raw-SAM conditional skewness or a dedicated mixed third moment, not inferred from covariance alone."
    )
    report.heading("Reproducibility")
    report.paragraph(
        "The report uses output/arm_stats.nc from the current PDF-10 full-LES-advance ARM run, output/arm_les_advance.nc for the supplied state, and the shared ARM 3-D SAM recreation.  Raw elapsed time equals CLUBB absolute time minus 41,400 s.  The JSON data artifact records every selected plane and component correlation.")
    report.code_file(Path(__file__), language="python", caption="Reproduction script (also copied into this bundle)")
    (report.stage / "data").mkdir(exist_ok=True)
    shutil.copy2(data_path, report.stage / "data" / data_path.name)
    report.publish()


def main() -> int:
    if not STATS_PATH.is_file() or not LES_PATH.is_file():
        raise FileNotFoundError("Expected current ARM PDF-10 full-LES-advance output is missing under output/.")
    with tempfile.TemporaryDirectory(prefix="clubb_pdf10_rt_thl_report_") as temporary:
        work_dir = Path(temporary)
        stats = xr.open_dataset(STATS_PATH, decode_times=False)
        les = xr.open_dataset(LES_PATH, decode_times=False)
        try:
            selection = _render_rt_thl_mosaic(work_dir, stats)
            projections = _render_projection_comparison(work_dir, stats)
            context = _render_forced_moment_context(work_dir, les)
            _publish(work_dir, selection, projections, context)
        finally:
            stats.close()
            les.close()
    print(f"Published doc/reports/{REPORT_ID}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
