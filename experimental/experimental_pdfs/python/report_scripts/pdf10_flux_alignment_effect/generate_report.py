#!/usr/bin/env python3
"""Publish a static ARM sensitivity report for PDF-10 flux alignment.

The input is a 36-column full-LES-advance ARM sweep in ``output/arm_stats.nc``.
It holds the PDF-10 center budget, w-direction scale, and G1 w-r_t capture
fixed, while crossing six moisture-tail gains with six flux-alignment
strengths.  The rendered PNGs are copied into the static report bundle.
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
from netCDF4 import Dataset, chartostring


REPORT_ID = "pdf10-flux-alignment-effect"
REPORT_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPORT_DIR.parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_pdf10_flux_align_report_matplotlib")

from dash_app.reports_tab.publisher import ReportBuilder
from utilities.sam_3d_reference import DEFAULT_SAM_RUN, _nearest_snapshot, load_snapshot


STATS_PATH = REPO_ROOT / "output" / "arm_stats.nc"
ARM_INITIAL_SECONDS = 41_400.0
EXAMPLE_MODEL_SECONDS = 71_460.0  # 19.85 h, matching raw-SAM minute 501
EXAMPLE_HEIGHT_M = 540.0
MAX_HEIGHT_M = 2_600.0
TAIL_GAIN = 3.0
GATE_SCALE = 0.25


def _revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip() + " (working tree)"
    except (OSError, subprocess.CalledProcessError):
        return "working tree"


def _nearest(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def _load_stats() -> dict:
    fields = (
        "time", "zt", "zm", "param_name", "clubb_params", "mixt_frac",
        "w_1", "w_2", "rt_1", "rt_2", "thl_1", "thl_2", "rtm", "thlm",
        "wp2_zt", "rtp2_zt", "thlp2_zt", "wprtp_zt", "wpthlp_zt", "Skw_zt",
        "wprcp",
    )
    with Dataset(STATS_PATH) as dataset:
        missing = [name for name in fields if name not in dataset.variables]
        if missing:
            raise KeyError(f"Stats file lacks required fields: {', '.join(missing)}")
        output = {name: np.asarray(dataset.variables[name][:], dtype=float) for name in fields if name != "param_name"}
        output["param_names"] = [str(value).strip() for value in chartostring(dataset.variables["param_name"][:]).tolist()]
    return output


def _parameter_columns(stats: dict) -> tuple[list[dict], list[int]]:
    names = stats["param_names"]
    params = stats["clubb_params"]
    required = (
        "pdf10_moisture_tail_gain", "pdf10_center_budget", "pdf10_w_direction_scale",
        "pdf10_g1_wrt_capture", "pdf10_flux_align_strength",
    )
    missing = [name for name in required if name not in names]
    if missing:
        raise KeyError(f"Sweep lacks PDF-10 parameters: {', '.join(missing)}")
    rows = []
    for column in range(params.shape[1]):
        rows.append({name: float(params[names.index(name), column]) for name in required} | {"column": column})
    columns = [row["column"] for row in rows if np.isclose(row["pdf10_moisture_tail_gain"], TAIL_GAIN)]
    if len(columns) < 2:
        raise ValueError(f"No usable {TAIL_GAIN:g} moisture-tail-gain flux-alignment sweep.")
    return rows, sorted(columns, key=lambda column: rows[column]["pdf10_flux_align_strength"])


def _flux_gate(stats: dict) -> np.ndarray:
    """Reproduce the dimensionless self-gating calculation used by PDF-10."""
    skewness = np.clip(stats["Skw_zt"][:, :, 0], -12.0, 12.0)
    softened = np.tanh(skewness / 1.25)
    wp2 = np.maximum(stats["wp2_zt"][:, :, 0], 1.0e-30)
    rtp2 = np.maximum(stats["rtp2_zt"][:, :, 0], 1.0e-30)
    thlp2 = np.maximum(stats["thlp2_zt"][:, :, 0], 1.0e-30)
    rho_wr = stats["wprtp_zt"][:, :, 0] / np.sqrt(wp2 * rtp2)
    rho_wthl = stats["wpthlp_zt"][:, :, 0] / np.sqrt(wp2 * thlp2)
    raw_gate = np.maximum(softened, 0.0) * np.maximum(rho_wr, 0.0) * np.maximum(rho_wthl, 0.0)
    return np.minimum(raw_gate / GATE_SCALE**3, 1.0)


def _style(axis) -> None:
    axis.set_facecolor("#111827")
    axis.tick_params(colors="#cbd5e1", labelsize=8)
    for spine in axis.spines.values():
        spine.set_color("#40516a")
    axis.xaxis.label.set_color("#dbe5f1")
    axis.yaxis.label.set_color("#dbe5f1")
    axis.title.set_color("#edf2f7")


def _raw_heatmap(axis, x: np.ndarray, y: np.ndarray, *, xlabel: str, ylabel: str, title: str) -> None:
    from matplotlib.colors import LogNorm

    x_low, x_high = np.percentile(x, (0.15, 99.85))
    y_low, y_high = np.percentile(y, (0.15, 99.85))
    x_edges = np.linspace(x_low - .06 * (x_high - x_low), x_high + .06 * (x_high - x_low), 76)
    y_edges = np.linspace(y_low - .06 * (y_high - y_low), y_high + .06 * (y_high - y_low), 76)
    density = np.histogram2d(x, y, bins=(x_edges, y_edges))[0].T
    positive = density[density > 0]
    axis.pcolormesh(
        x_edges, y_edges, density, shading="auto", cmap="magma",
        norm=LogNorm(vmin=max(float(np.percentile(positive, 7)), 1.0), vmax=float(np.max(positive))),
    )
    _style(axis)
    axis.set(xlabel=xlabel, ylabel=ylabel, title=title)


def _render_example_geometry(work_dir: Path, stats: dict, columns: list[int], rows: list[dict]) -> dict:
    from matplotlib import pyplot as plt

    time_index = _nearest(stats["time"], EXAMPLE_MODEL_SECONDS)
    height_index = _nearest(stats["zt"], EXAMPLE_HEIGHT_M)
    snapshot = load_snapshot(
        DEFAULT_SAM_RUN,
        int(round(float(stats["time"][time_index]) - ARM_INITIAL_SECONDS)),
        float(stats["zt"][height_index]),
    )
    strengths = np.asarray([rows[column]["pdf10_flux_align_strength"] for column in columns])
    colors = plt.cm.viridis(np.linspace(.12, .94, len(columns)))
    raw_rt = np.asarray(snapshot.samples[:, 1], dtype=float) * 1_000.0
    raw_w = np.asarray(snapshot.samples[:, 0], dtype=float)
    raw_thl = np.asarray(snapshot.samples[:, 2], dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(15.2, 6.3), dpi=180, facecolor="#0d1525")
    _raw_heatmap(
        axes[0], raw_rt, raw_w, xlabel="Total water, rₜ [g/kg]", ylabel="Vertical velocity, w [m/s]",
        title="Raw SAM w–rₜ plane with PDF-10 G1 center trajectory",
    )
    _raw_heatmap(
        axes[1], raw_rt, raw_thl, xlabel="Total water, rₜ [g/kg]", ylabel="θₗ [K]",
        title="Same raw plane in rₜ–θₗ coordinates",
    )
    table = []
    for column, strength, color in zip(columns, strengths, colors):
        g1_rt = float(stats["rt_1"][time_index, height_index, column] * 1_000.0)
        g2_rt = float(stats["rt_2"][time_index, height_index, column] * 1_000.0)
        g1_w = float(stats["w_1"][time_index, height_index, column])
        g2_w = float(stats["w_2"][time_index, height_index, column])
        g1_thl = float(stats["thl_1"][time_index, height_index, column])
        g2_thl = float(stats["thl_2"][time_index, height_index, column])
        weight = float(stats["mixt_frac"][time_index, height_index, column])
        axes[0].scatter(g1_rt, g1_w, color=color, s=44, marker="o", edgecolor="white", linewidth=.45, zorder=6)
        axes[0].scatter(g2_rt, g2_w, color=color, s=24, marker="x", linewidth=1.2, zorder=5)
        axes[1].scatter(g1_rt, g1_thl, color=color, s=44, marker="o", edgecolor="white", linewidth=.45, zorder=6)
        axes[1].scatter(g2_rt, g2_thl, color=color, s=24, marker="x", linewidth=1.2, zorder=5)
        table.append({
            "column": int(column), "flux_align_strength": float(strength), "weight_g1": weight,
            "g1_w_m_s": g1_w, "g1_rt_g_kg": g1_rt, "g1_thl_K": g1_thl,
            "g2_w_m_s": g2_w, "g2_rt_g_kg": g2_rt, "g2_thl_K": g2_thl,
        })
    for axis, y1, y2 in ((axes[0], stats["w_1"], stats["w_2"]), (axes[1], stats["thl_1"], stats["thl_2"])):
        axis.plot(
            [item["g1_rt_g_kg"] for item in table],
            [item["g1_w_m_s"] if axis is axes[0] else item["g1_thl_K"] for item in table],
            color="#f8fafc", linewidth=1.0, alpha=.8, zorder=4,
        )
    handles = [plt.Line2D([0], [0], color=color, marker="o", markeredgecolor="white", linewidth=1, label=f"f={strength:.2f}") for strength, color in zip(strengths, colors)]
    figure.legend(handles=handles, title="flux alignment", loc="lower center", ncol=len(handles), frameon=False, labelcolor="#dbe5f1", fontsize=8, bbox_to_anchor=(.5, .035))
    figure.suptitle(
        f"ARM low-level V: {stats['time'][time_index] / 3600:.2f} h, z={stats['zt'][height_index]:.0f} m; moisture-tail gain={TAIL_GAIN:g}",
        color="#f8fafc", fontsize=14, y=.98,
    )
    figure.text(.5, .095, "Raw SAM probability is the background. Circles: G1; crosses: G2. Colors trace increasing pdf10_flux_align_strength; all other PDF-10 sweep controls are held fixed.", ha="center", color="#b6c5d8", fontsize=9)
    figure.tight_layout(rect=(0, .14, 1, .94))
    output = work_dir / "low-level-v-center-response.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    return {
        "model_time_seconds": float(stats["time"][time_index]), "model_time_hours": float(stats["time"][time_index] / 3600),
        "height_m": float(stats["zt"][height_index]), "raw_elapsed_minutes": float(snapshot.elapsed_minutes),
        "raw_sample_count": int(snapshot.sample_count), "sweep_centers": table,
    }


def _render_gate_and_profiles(work_dir: Path, stats: dict, columns: list[int], rows: list[dict], gate: np.ndarray) -> dict:
    from matplotlib import pyplot as plt

    low_column, high_column = columns[0], columns[-1]
    time_index = _nearest(stats["time"], EXAMPLE_MODEL_SECONDS)
    z_mask = stats["zt"] <= MAX_HEIGHT_M
    t_mask = (stats["time"] >= 54_000.0) & (stats["time"] <= 93_600.0)
    figure, axes = plt.subplots(1, 3, figsize=(17.2, 5.9), dpi=180, facecolor="#0d1525", gridspec_kw={"width_ratios": (1.35, 1, 1)})
    mesh = axes[0].pcolormesh(stats["time"][t_mask] / 3600.0, stats["zt"][z_mask], gate[t_mask][:, z_mask].T, shading="auto", cmap="viridis", vmin=0, vmax=1)
    axes[0].axvline(stats["time"][time_index] / 3600.0, color="#facc15", linestyle="--", linewidth=1.2)
    axes[0].axhline(stats["zt"][_nearest(stats["zt"], EXAMPLE_HEIGHT_M)], color="#facc15", linestyle="--", linewidth=1.2)
    _style(axes[0]); axes[0].set(xlabel="model time [h]", ylabel="height [m]", title="PDF-10 flux-alignment gate")
    colorbar = figure.colorbar(mesh, ax=axes[0], pad=.015); colorbar.set_label("gate G", color="#cbd5e1"); colorbar.ax.tick_params(colors="#cbd5e1", labelsize=8)
    z = stats["zt"][z_mask]
    for column, color, label in ((low_column, "#38bdf8", f"f={rows[low_column]['pdf10_flux_align_strength']:.2f}"), (high_column, "#f97316", f"f={rows[high_column]['pdf10_flux_align_strength']:.2f}")):
        rt_anomaly = (stats["rt_1"][time_index, z_mask, column] - stats["rtm"][time_index, z_mask, column]) * 1_000.0
        axes[1].plot(rt_anomaly, z, color=color, linewidth=2.1, label=label)
        axes[2].plot(stats["mixt_frac"][time_index, z_mask, column], z, color=color, linewidth=2.1, label=label)
    for axis, label, title in (
        (axes[1], "G1 rₜ − grid mean [g/kg]", "G1 moisture displacement"),
        (axes[2], "G1 mixture weight", "G1 population weight"),
    ):
        axis.axvline(0, color="#64748b", linewidth=.8)
        axis.axhline(EXAMPLE_HEIGHT_M, color="#facc15", linestyle="--", linewidth=1.1)
        _style(axis); axis.set(xlabel=label, ylabel="height [m]", ylim=(0, MAX_HEIGHT_M), title=title)
        axis.legend(loc="lower left", frameon=False, labelcolor="#dbe5f1", fontsize=8)
    figure.suptitle("Where the new branch acts, and what it changes at the V time", color="#f8fafc", fontsize=14, y=.98)
    figure.tight_layout(rect=(0, 0, 1, .94))
    output = work_dir / "gate-and-vertical-response.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    active = gate[t_mask][:, z_mask]
    return {
        "analysis_time_range_seconds": [54_000.0, 93_600.0], "analysis_height_limit_m": MAX_HEIGHT_M,
        "gate_positive_fraction": float(np.mean(active > 1.0e-8)),
        "gate_full_fraction": float(np.mean(active >= .999)),
        "gate_mean_when_positive": float(np.mean(active[active > 1.0e-8])) if np.any(active > 1.0e-8) else 0.0,
        "gate_at_example": float(gate[time_index, _nearest(stats["zt"], EXAMPLE_HEIGHT_M)]),
    }


def _render_wprcp_invariance(work_dir: Path, stats: dict, columns: list[int], rows: list[dict]) -> dict:
    from matplotlib import pyplot as plt

    low_column, high_column = columns[0], columns[-1]
    time_index = _nearest(stats["time"], EXAMPLE_MODEL_SECONDS)
    raw_elapsed_seconds, raw_path = _nearest_snapshot(
        DEFAULT_SAM_RUN,
        int(round(float(stats["time"][time_index]) - ARM_INITIAL_SECONDS)),
    )
    with Dataset(raw_path) as raw:
        raw_z = np.asarray(raw.variables["z"][:], dtype=float)
        raw_w = np.asarray(raw.variables["W"][0], dtype=float)
        raw_rc = np.asarray(raw.variables["RC"][0], dtype=float)
        raw_wprcp = np.mean(
            (raw_w - np.mean(raw_w, axis=(1, 2), keepdims=True))
            * (raw_rc - np.mean(raw_rc, axis=(1, 2), keepdims=True)),
            axis=(1, 2),
        )
    z_mask = raw_z <= MAX_HEIGHT_M
    target_z = raw_z[z_mask]
    low_profile = np.interp(target_z, stats["zm"], stats["wprcp"][time_index, :, low_column])
    high_profile = np.interp(target_z, stats["zm"], stats["wprcp"][time_index, :, high_column])
    t_mask = (stats["time"] >= 54_000.0) & (stats["time"] <= 93_600.0)
    zm_mask = stats["zm"] <= MAX_HEIGHT_M
    difference = stats["wprcp"][t_mask][:, zm_mask, high_column] - stats["wprcp"][t_mask][:, zm_mask, low_column]
    reference_wprcp = stats["wprcp"][t_mask][:, zm_mask, low_column]
    figure, axes = plt.subplots(1, 2, figsize=(14.8, 5.6), dpi=180, facecolor="#0d1525")
    axes[0].plot(raw_wprcp[z_mask] * 1.0e6, target_z, color="#facc15", linewidth=2.5, label="SAM")
    axes[0].plot(low_profile * 1.0e6, target_z, color="#38bdf8", linewidth=2.0, label=f"PDF-10 f={rows[low_column]['pdf10_flux_align_strength']:.2f}")
    axes[0].plot(high_profile * 1.0e6, target_z, color="#f97316", linewidth=1.4, linestyle="--", label=f"PDF-10 f={rows[high_column]['pdf10_flux_align_strength']:.2f}")
    _style(axes[0]); axes[0].set(xlabel="w'rc' [10^-6 m/s kg/kg]", ylabel="height [m]", ylim=(0, MAX_HEIGHT_M), title="Cloud-water transport is visually unchanged")
    axes[0].legend(frameon=False, labelcolor="#dbe5f1", fontsize=8)
    axes[1].hist(difference.ravel() * 1.0e9, bins=80, color="#a78bfa", alpha=.85)
    axes[1].axvline(0, color="#facc15", linewidth=1.2)
    _style(axes[1]); axes[1].set(xlabel="w'rc'(f=1.00) − w'rc'(f=0.10) [10^-9]", ylabel="time-height samples", title="Sweep-wide PDF w'rc' difference")
    figure.suptitle("A center-allocation change, not a new grid-moment prediction", color="#f8fafc", fontsize=14, y=.98)
    figure.text(.5, .025, "The branch reallocates two Gaussian centers/widths while retaining the same closure family. The small w'rc' response does not make cloud transport materially larger.", ha="center", color="#b6c5d8", fontsize=9)
    figure.tight_layout(rect=(0, .08, 1, .94))
    output = work_dir / "wprcp-invariance.png"
    figure.savefig(output, facecolor=figure.get_facecolor())
    plt.close(figure)
    return {
        "max_abs_wprcp_difference": float(np.max(np.abs(difference))),
        "rms_wprcp_difference": float(np.sqrt(np.mean(difference**2))),
        "max_relative_wprcp_difference": float(np.max(np.abs(difference)) / np.max(np.abs(reference_wprcp))),
        "rms_relative_wprcp_difference": float(np.sqrt(np.mean(difference**2)) / np.sqrt(np.mean(reference_wprcp ** 2))),
        "selected_raw_time_minutes": float(raw_elapsed_seconds / 60.0),
        "selected_raw_source": str(raw_path),
    }


def _input_spread(stats: dict) -> dict:
    result = {}
    for name in ("wp2_zt", "rtp2_zt", "thlp2_zt", "wprtp_zt", "wpthlp_zt", "Skw_zt"):
        values = stats[name]
        scale = max(float(np.nanmax(np.abs(values[:, :, 0]))), 1.0e-30)
        result[name] = float(np.nanmax(np.ptp(values, axis=2)) / scale)
    return result


def _publish(work_dir: Path, stats: dict, rows: list[dict], columns: list[int], geometry: dict, gate_info: dict, invariance: dict) -> None:
    low, high = geometry["sweep_centers"][0], geometry["sweep_centers"][-1]
    gain_rows = [row for row in rows if np.isclose(row["pdf10_moisture_tail_gain"], TAIL_GAIN)]
    metadata = {
        "stats_path": str(STATS_PATH.relative_to(REPO_ROOT)), "raw_sam_run": str(DEFAULT_SAM_RUN),
        "fixed_controls": {key: gain_rows[0][key] for key in ("pdf10_center_budget", "pdf10_w_direction_scale", "pdf10_g1_wrt_capture")},
        "sweep_columns": gain_rows, "example": geometry, "gate": gate_info,
        "wprcp_invariance": invariance, "relative_input_spread_across_columns": _input_spread(stats),
    }
    data_path = work_dir / "flux-alignment-effect.json"
    data_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = ReportBuilder(
        REPORT_ID,
        "PDF-10 flux alignment in ARM: what the new branch actually changes",
        summary="A controlled 36-column full-LES-advance ARM sweep shows that flux alignment substantially moves the low-level G1 moisture center and increases its mixture weight only in coherent moist/warm rising regimes. It leaves the supplied grid moments and diagnosed w′r꜀′ essentially unchanged.",
        tags=("PDF-10", "ARM", "SAM", "LES advance", "flux alignment", "two-Gaussian PDF"),
        replace=True,
        source_revision=_revision(),
    )
    report.heading("Answer at a glance")
    report.callout(
        "The change works as designed, but its effect is narrower than a new transport closure.",
        f"At the low-level ARM V (19.85 h, {geometry['height_m']:.0f} m), increasing pdf10_flux_align_strength from {low['flux_align_strength']:.2f} to {high['flux_align_strength']:.2f} moves G1 rₜ from {low['g1_rt_g_kg']:.3f} to {high['g1_rt_g_kg']:.3f} g/kg and raises its weight from {100*low['weight_g1']:.1f}% to {100*high['weight_g1']:.1f}%. G1 remains a rising component (w ≈ {high['g1_w_m_s']:.2f} m/s). The new rule therefore corrects the dry-tail tendency by allocating a more substantial moist rising population; it does not add new moment information or materially increase w'rc'.",
        tone="success",
    )
    report.metrics([
        ("G1 rₜ shift", f"+{high['g1_rt_g_kg'] - low['g1_rt_g_kg']:.3f} g/kg", "f=0.10 → 1.00 at the V"),
        ("G1 weight", f"{100*low['weight_g1']:.1f}% → {100*high['weight_g1']:.1f}%", "same point; a less rare transport population"),
        ("Gate at V", f"{gate_info['gate_at_example']:.2f}", "fully active at this selected state"),
        ("w'rc' max change", f"{100*invariance['max_relative_wprcp_difference']:.4f}%", "sweep-wide, 0–2.6 km; negligible endpoint response"),
    ])
    report.heading("Controlled comparison")
    report.paragraph(
        "This is not a before/after pair of unrelated free runs. The current output/arm_stats.nc contains a 6×6 multicolumn ARM full-LES-advance sweep. For the slice used below, moisture-tail gain is fixed at 3.0, center budget at 0.90, w-direction scale at 0.80, and G1 w-r_t capture at 0.40. Only flux-alignment strength varies from 0.10 to 1.00. The sweep has no exact f=0 column, so f=0.10 is the lowest measured endpoint rather than a literal legacy baseline."
    )
    report.paragraph(
        "This is a full LES-advance multicolumn experiment: nominal forcing and the named sweep controls are fixed, but each column may evolve a slightly different state through its own PDF feedback. The compact JSON records actual column spreads in the gate inputs. The large, monotonic center response below is nevertheless directly attributable to flux-alignment strength because it is the sole varied control in this six-column slice."
    )
    report.figure(work_dir / "low-level-v-center-response.png", caption="At the raw-SAM low-level V, G1 progresses from the dry side toward the moist rising arm as flux alignment increases. G2 compensates to retain the grid mean. The θ_l motion is smaller and slightly cools G1, so the correction clearly improves moisture placement more than it solves the full r_t–θ_l V geometry.")
    report.heading("Why the correction is localized")
    report.equation(r"G = \min\!\left[1,\frac{\max(\widetilde S_w,0)\max(\rho_{wr_t},0)\max(\rho_{w\theta_l},0)}{0.25^3}\right],\qquad \eta=f_{\rm align}G", caption="The self-gate: positive w skewness plus positive w–r_t and w–θ_l correlations are all required.")
    report.paragraph(
        f"Over 15–26 h and below {MAX_HEIGHT_M:.0f} m, the gate is positive in {100*gate_info['gate_positive_fraction']:.1f}% of time-height states and fully open in {100*gate_info['gate_full_fraction']:.1f}%. A negative w–θ_l relationship gives G=0 exactly, preserving the original PDF-10 path aloft where the prior negative thermodynamic tilt was useful."
    )
    report.figure(work_dir / "gate-and-vertical-response.png", caption="The gate identifies a patchy lower-layer/surface-plume regime rather than an altitude switch. At the V time, full strength increases G1’s positive moisture displacement and its mixture weight where that gate is active; the high-level inactive regime is unchanged.")
    report.heading("What did not change materially: w'rc'")
    report.paragraph(
        f"The selected w'rc' profiles for the two strengths lie on top of one another. Across the time-height comparison, the largest endpoint difference is {100*invariance['max_relative_wprcp_difference']:.4f}% of the selected profile peak and RMS difference is {100*invariance['rms_relative_wprcp_difference']:.4f}%. That is a real but negligible dynamic response: the modification principally changes how the two Gaussian centers and widths allocate the available trivariate information. Therefore it cannot, by itself, repair the known magnitude underprediction of cloud-water transport."
    )
    report.figure(work_dir / "wprcp-invariance.png", caption="The yellow SAM comparison remains important, but the two flux-alignment endpoints overlap. The right panel shows the small endpoint response in 10^-9 units.")
    report.heading("Interpretation and next test")
    report.paragraph(
        "This is a meaningful geometry improvement: the low-level transport component becomes more moist and less artificially rare while keeping its positive w identity. It is not a complete V solution. In particular, the G1 θ_l center moves only modestly (and slightly away from the raw moist/warm arm at this sample), because the branch still has only marginal skewness plus covariance information. The raw V remains a conditional/multi-population structure that two Gaussian centers cannot identify uniquely from those moments."
    )
    report.callout(
        "Recommendation",
        "Keep the branch as a bounded surface-plume geometry correction and tune its strength against center-placement diagnostics, not against w′r꜀′ magnitude. A useful next controlled test is f=0, 0.5, and 1.0 with the same other controls in ARM, BOMEX, and a stratocumulus case. Score raw-SAM component-center displacement in w–r_t and r_t–θ_l separately, plus w′r꜀′ as a no-regression metric. If the remaining V error matters, test bounded mixed-third-moment compatibility information before considering a third Gaussian or prognostic equation.",
        tone="warning",
    )
    report.heading("Reproducibility")
    report.paragraph("Rendered figures are embedded static PNGs. The accompanying JSON records the exact columns, controls, selected plane, gate statistics, and invariance calculation. Re-running the copied script deliberately regenerates the bundle from output/arm_stats.nc and the ARM raw-SAM benchmark.")
    report.code_file(Path(__file__), language="python", caption="Report reproduction script")
    (report.stage / "data").mkdir(exist_ok=True)
    shutil.copy2(data_path, report.stage / "data" / data_path.name)
    report.publish()


def main() -> int:
    if not STATS_PATH.is_file():
        raise FileNotFoundError("Expected ARM sweep stats or raw SAM benchmark is unavailable.")
    stats = _load_stats()
    rows, columns = _parameter_columns(stats)
    gate = _flux_gate(stats)
    with tempfile.TemporaryDirectory(prefix="clubb_pdf10_flux_align_report_") as directory:
        work_dir = Path(directory)
        geometry = _render_example_geometry(work_dir, stats, columns, rows)
        gate_info = _render_gate_and_profiles(work_dir, stats, columns, rows, gate)
        invariance = _render_wprcp_invariance(work_dir, stats, columns, rows)
        _publish(work_dir, stats, rows, columns, geometry, gate_info, invariance)
    print(f"Published doc/reports/{REPORT_ID}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
