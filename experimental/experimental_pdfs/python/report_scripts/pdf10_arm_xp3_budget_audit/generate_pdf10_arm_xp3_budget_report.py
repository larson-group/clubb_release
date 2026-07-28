#!/usr/bin/env python3
"""Publish the PDF-10 ARM scalar-third-moment budget audit.

The report reads dedicated SCM artifacts rather than ``output/``.  It is
therefore stable when interactive Dash runs later replace the usual ARM files.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dash_app.reports_tab.publisher import ReportBuilder


ARTIFACT_ROOT = ROOT / "agent_artifacts" / "xp3-budget-arm-20260724"
RAW_SAM = (
    ROOT
    / "input/les_and_clubb_benchmark_runs/sam_benchmark_runs/JULY_2017/ARM_96x96x110"
    / "GCSSARM_96x96x110_67m_40m_1s.nc"
)
RUNS = {
    "before": ARTIFACT_ROOT / "les-lower-no-xp3" / "arm_stats.nc",
    "corrected": ARTIFACT_ROOT / "les-lower-no-xp3-zt-pdf-flux" / "arm_stats.nc",
}
REPORT_ID = "pdf10-arm-xp3-budget-audit"
REPORT_TITLE = "PDF-10 ARM scalar-third-moment budget: first resolved transport audit"
SAM_CLOCK_OFFSET_SECONDS = 11.5 * 3600.0
SNAPSHOT_SECONDS = 66_600.0
HEIGHT_LIMIT_M = 2_500.0


def _revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _field(ds: xr.Dataset, name: str) -> np.ndarray:
    value = ds[name]
    vertical = "zm" if "zm" in value.dims else "zt" if "zt" in value.dims else "z"
    for dimension in tuple(value.dims):
        if dimension not in {"time", vertical}:
            value = value.isel({dimension: 0})
    return np.asarray(value.values, dtype=float)


def _raw_field(ds: xr.Dataset, name: str) -> np.ndarray:
    value = ds[name]
    for dimension in tuple(value.dims):
        if dimension not in {"time", "z"}:
            value = value.isel({dimension: 0})
    return np.asarray(value.values, dtype=float)


def _style(axis: plt.Axes) -> None:
    axis.set_facecolor("#111827")
    axis.tick_params(colors="#cbd5e1", labelsize=8)
    axis.xaxis.label.set_color("#dbe5f1")
    axis.yaxis.label.set_color("#dbe5f1")
    axis.title.set_color("#f8fafc")
    axis.grid(color="#334155", linewidth=0.55, alpha=0.7)
    for spine in axis.spines.values():
        spine.set_color("#40516a")


def _save(figure: plt.Figure, path: Path) -> None:
    figure.tight_layout()
    figure.savefig(path, dpi=100, facecolor=figure.get_facecolor())
    plt.close(figure)


def _interp_to_clubb(raw_values: np.ndarray, raw_time: np.ndarray, raw_z: np.ndarray,
                     clubb_time: np.ndarray, clubb_z: np.ndarray) -> np.ndarray:
    temporal = np.empty((len(clubb_time), len(raw_z)), dtype=float)
    for k in range(len(raw_z)):
        temporal[:, k] = np.interp(clubb_time, raw_time, raw_values[:, k])
    output = np.empty((len(clubb_time), len(clubb_z)), dtype=float)
    for t in range(len(clubb_time)):
        output[t] = np.interp(clubb_z, raw_z, temporal[t], left=np.nan, right=np.nan)
    return output


def _load() -> tuple[dict[str, xr.Dataset], xr.Dataset]:
    missing = [str(path) for path in (*RUNS.values(), RAW_SAM) if not path.is_file()]
    if missing:
        raise FileNotFoundError("required audit artifact(s) missing: " + ", ".join(missing))
    return ({name: xr.open_dataset(path, decode_times=False) for name, path in RUNS.items()},
            xr.open_dataset(RAW_SAM, decode_times=False))


def _metrics(runs: dict[str, xr.Dataset], raw: xr.Dataset) -> dict[str, float]:
    corrected = runs["corrected"]
    time = np.asarray(corrected["time"].values, dtype=float)
    zt = np.asarray(corrected["zt"].values, dtype=float)
    raw_time = np.asarray(raw["time"].values, dtype=float) * 60.0 + SAM_CLOCK_OFFSET_SECONDS
    raw_z = np.asarray(raw["z"].values, dtype=float)
    ti = int(np.argmin(np.abs(time - SNAPSHOT_SECONDS)))
    mask = zt <= HEIGHT_LIMIT_M
    result: dict[str, float] = {"snapshot_seconds": float(time[ti])}
    for variable, raw_name in (("rtp3", "RTP3"), ("thlp3", "THLP3")):
        reference = _interp_to_clubb(_raw_field(raw, raw_name), raw_time, raw_z, time, zt)
        valid_reference = mask & np.isfinite(reference[ti])
        result[f"sam_{variable}_snapshot_abs_peak"] = float(
            np.max(np.abs(reference[ti, valid_reference]))
        )
        for label, ds in runs.items():
            value = _field(ds, variable)
            valid = mask & np.isfinite(value[ti]) & np.isfinite(reference[ti])
            result[f"{label}_{variable}_snapshot_abs_peak"] = float(
                np.max(np.abs(value[ti, valid]))
            )
            result[f"{label}_{variable}_snapshot_rmse"] = float(
                np.sqrt(np.mean((value[ti, valid] - reference[ti, valid]) ** 2))
            )
            result[f"{label}_{variable}_snapshot_correlation"] = float(
                np.corrcoef(value[ti, valid], reference[ti, valid])[0, 1]
            )
    for variable in ("rtp3", "thlp3"):
        budget_names = [f"{variable}_{suffix}" for suffix in ("ta", "tp", "ac", "df", "dp")]
        rms = {
            name: float(np.sqrt(np.nanmean(_field(corrected, name)[:, mask] ** 2)))
            for name in budget_names
        }
        result.update({f"{name}_rms": value for name, value in rms.items()})
        result[f"{variable}_residual_ratio"] = float(
            np.sqrt(np.nanmean(_field(corrected, f"{variable}_rs")[:, mask] ** 2))
            / max(sum(rms.values()), np.finfo(float).tiny)
        )
    return result


def _plot_profiles(work: Path, runs: dict[str, xr.Dataset], raw: xr.Dataset) -> None:
    corrected = runs["corrected"]
    time = np.asarray(corrected["time"].values, dtype=float)
    zt = np.asarray(corrected["zt"].values, dtype=float)
    raw_time = np.asarray(raw["time"].values, dtype=float) * 60.0 + SAM_CLOCK_OFFSET_SECONDS
    raw_z = np.asarray(raw["z"].values, dtype=float)
    ti = int(np.argmin(np.abs(time - SNAPSHOT_SECONDS)))
    ri = int(np.argmin(np.abs(raw_time - time[ti])))
    mask = zt <= HEIGHT_LIMIT_M
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 6.4), sharey=True, facecolor="#0d1525")
    for axis, variable, raw_name, unit in (
        (axes[0], "rtp3", "RTP3", r"r_t'^3$ [(kg kg$^{-1}$)$^3$]"),
        (axes[1], "thlp3", "THLP3", r"theta_l'^3$ [K$^3$]"),
    ):
        axis.plot(_raw_field(raw, raw_name)[ri], raw_z, color="#facc15", linewidth=2.2, label="raw SAM")
        axis.plot(_field(runs["before"], variable)[ti, mask], zt[mask], color="#64748b",
                  linewidth=1.5, linestyle="--", label="before flux wiring")
        axis.plot(_field(corrected, variable)[ti, mask], zt[mask], color="#38bdf8",
                  linewidth=1.8, label="corrected PDF flux")
        axis.axvline(0.0, color="#94a3b8", linewidth=0.75)
        axis.set(xlabel=unit, ylim=(0.0, HEIGHT_LIMIT_M), title=variable)
        _style(axis)
    axes[0].set_ylabel("height [m]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, labelcolor="#e2e8f0", fontsize=8)
    fig.suptitle(f"ARM at {time[ti] / 3600:.2f} h: SAM-clock aligned scalar third moments", color="#f8fafc", fontsize=14)
    fig.subplots_adjust(bottom=0.15, top=0.90, wspace=0.24)
    fig.savefig(work / "profiles-before-after.png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)


def _plot_budget_profiles(work: Path, corrected: xr.Dataset) -> None:
    time = np.asarray(corrected["time"].values, dtype=float)
    zt = np.asarray(corrected["zt"].values, dtype=float)
    ti = int(np.argmin(np.abs(time - SNAPSHOT_SECONDS)))
    mask = zt <= HEIGHT_LIMIT_M
    colors = {"bt": "#f8fafc", "ta": "#38bdf8", "tp": "#facc15", "ac": "#f97316", "df": "#a78bfa", "dp": "#ec4899"}
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 6.4), sharey=True, facecolor="#0d1525")
    for axis, variable, unit in (
        (axes[0], "rtp3", r"tendency [(kg kg$^{-1}$)$^3$ s$^{-1}$]"),
        (axes[1], "thlp3", r"tendency [K$^3$ s$^{-1}$]"),
    ):
        for suffix in ("bt", "ta", "tp", "ac", "df", "dp"):
            axis.plot(_field(corrected, f"{variable}_{suffix}")[ti, mask], zt[mask],
                      color=colors[suffix], linewidth=1.45 if suffix != "bt" else 2.0, label=suffix)
        axis.axvline(0.0, color="#94a3b8", linewidth=0.75)
        axis.set(xlabel=unit, ylim=(0.0, HEIGHT_LIMIT_M), title=f"{variable} budget")
        _style(axis)
    axes[0].set_ylabel("height [m]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, labelcolor="#e2e8f0", fontsize=8)
    fig.suptitle(f"Corrected PDF-10 tendency balance at {time[ti] / 3600:.2f} h", color="#f8fafc", fontsize=14)
    fig.subplots_adjust(bottom=0.15, top=0.90, wspace=0.24)
    fig.savefig(work / "budget-profiles.png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)


def _plot_time_height(work: Path, corrected: xr.Dataset) -> None:
    time_hours = np.asarray(corrected["time"].values, dtype=float) / 3600.0
    zt = np.asarray(corrected["zt"].values, dtype=float)
    mask = zt <= HEIGHT_LIMIT_M
    panels = (("rtp3_ta", "PDF transport: r_t'^3"), ("rtp3_dp", "tau damping: r_t'^3"),
              ("thlp3_ta", "PDF transport: theta_l'^3"), ("thlp3_dp", "tau damping: theta_l'^3"))
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.5), sharex=True, sharey=True, facecolor="#0d1525")
    for axis, (name, title) in zip(axes.flat, panels):
        values = _field(corrected, name)[:, mask]
        limit = float(np.nanpercentile(np.abs(values), 99.0))
        image = axis.pcolormesh(time_hours, zt[mask], values.T, shading="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
        axis.set(title=title, ylim=(0.0, HEIGHT_LIMIT_M))
        _style(axis)
        colorbar = fig.colorbar(image, ax=axis, pad=0.015)
        colorbar.ax.tick_params(labelsize=7, colors="#cbd5e1")
    axes[0, 0].set_ylabel("height [m]")
    axes[1, 0].set_ylabel("height [m]")
    axes[1, 0].set_xlabel("time [h]")
    axes[1, 1].set_xlabel("time [h]")
    fig.suptitle("Resolved PDF transport is spatially organized, but it is not yet a complete source closure", color="#f8fafc", fontsize=14)
    fig.subplots_adjust(top=0.91, wspace=0.22, hspace=0.2)
    fig.savefig(work / "transport-and-damping-time-height.png", dpi=100, facecolor=fig.get_facecolor())
    plt.close(fig)


def _publish(work: Path, metrics: dict[str, float]) -> None:
    payload_path = work / "budget-audit-data.json"
    payload_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = ReportBuilder(
        REPORT_ID,
        REPORT_TITLE,
        summary="A corrected PDF-10 turbulent-third-moment flux now enters the ARM rtp3/thlp3 predictor, and the printed budget closes to roundoff.  The new transport restores substantial amplitude but does not yet place the scalar-third-moment structures correctly, pointing next to missing forcing covariance and PDF-geometry fidelity rather than numerical residuals.",
        tags=("PDF-10", "ARM", "rtp3", "thlp3", "budget", "LES advance"),
        replace=True,
        source_revision=_revision(),
    )
    report.heading("Answer at a glance")
    report.callout(
        "The budget is now scientifically usable, but the closure is not complete.",
        "The initial transport diagnostic was accidentally zero because it read the optional momentum-level PDF parameters.  PDF-10 normally uses one closure call, so that container is unset.  The corrected implementation evaluates <w'x'^3> from the thermodynamic-level PDF state and maps it to the flux grid.  Its discrete residual is roundoff-sized, and the turbulent-transport tendency is now an active leading term.",
        tone="success",
    )
    report.metrics([
        ("r_t residual", f"{metrics['rtp3_residual_ratio']:.1e}", "RMS residual / sum of RMS named tendencies"),
        ("theta_l residual", f"{metrics['thlp3_residual_ratio']:.1e}", "same discrete-closure check"),
        ("r_t transport RMS", f"{metrics['rtp3_ta_rms']:.2e}", "new conservative PDF transport tendency"),
        ("theta_l transport RMS", f"{metrics['thlp3_ta_rms']:.2e}", "new conservative PDF transport tendency"),
    ])
    report.heading("Configuration and comparison discipline")
    report.paragraph(
        "Both retained ARM SCM runs use PDF-10 and LES advance for means, scalar second moments/covariance, velocity moments, and horizontal winds, but deliberately do not advance rtp3 or thlp3.  Thus the scalar-third-moment predictor is tested against SAM while lower moments are supplied.  Raw SAM timestamps are aligned to CLUBB seconds since midnight by adding the 11:30 UTC (41400 s) archive origin; earlier quick-look comparisons that treated SAM minutes as absolute seconds would have compared the wrong physical time."
    )
    report.code(
        "-les_advance wp23,xm,xp2,windm\n"
        "iiPDF_type=10, l_predict_upwp_vpwp=false\n"
        "pdf10_moisture_tail_gain=2.51273, pdf10_center_budget=0.95\n"
        "pdf10_w_direction_scale=0.68, pdf10_g1_wrt_capture=0.131891\n"
        "pdf10_plume_structure_strength=0.168881, c_K3=0.025, nu3=1.0",
        language="text", caption="Dedicated retained SCM configuration",
    )
    report.heading("What changed in the equation")
    report.equation(
        r"\frac{d\overline{x'^3}}{dt}=T_{ta}+T_{tp}+T_{ac}+T_{ma}+T_{df}+T_{dp}",
        caption="Named terms now output for x = r_t and theta_l; the printed residual is the discrete difference between the stored tendency and their sum.",
    )
    report.paragraph(
        "For each Gaussian component, PDF-10 evaluates the exact frozen-state flux contribution d_w(d_x^3+3d_x V_x)+3C_wx(d_x^2+V_x).  The two component contributions are mixture weighted, evaluated on thermodynamic levels, interpolated to the momentum/flux grid, and inserted as the conservative density-weighted divergence.  Mean advection and diffusion remain in the existing thermodynamic-level tridiagonal solve; ARM has zero resolved mean vertical velocity here, so its named mean-advection term is correctly zero."
    )
    report.heading("Evidence: the transport fix changes amplitude, not yet placement")
    report.figure(work / "profiles-before-after.png", caption="At 18.50 h, the corrected transport gives r_t'^3 an absolute peak of %.2e, close to the raw-SAM peak of %.2e, compared with %.2e before the wiring fix.  However, profile correlation is still only %.2f.  theta_l'^3 similarly gains the negative excursion but remains spatially biased.  The result is a real source contribution, not a complete closure." % (
        metrics["corrected_rtp3_snapshot_abs_peak"], metrics["sam_rtp3_snapshot_abs_peak"],
        metrics["before_rtp3_snapshot_abs_peak"],
        metrics["corrected_rtp3_snapshot_correlation"],
    ))
    report.figure(work / "budget-profiles.png", caption="The corrected tendencies show where transport competes with production, accumulation, diffusion, and tau damping.  The r_t transport contribution is no longer identically zero; theta_l transport is active but not sufficient to reproduce the observed negative third-moment structure by itself.")
    report.figure(work / "transport-and-damping-time-height.png", caption="Time-height fields distinguish a physically structured transport tendency from a numerical closure error.  The residual remains negligible throughout; transport and damping are the prominent competing terms, while mean advection is zero in this anelastic ARM configuration.")
    report.heading("Interpretation")
    report.paragraph(
        "The correction validates the numerical architecture: the flux is on the correct staggered grid, its boundary flux is closed, and the discrete equation closes.  It also sharpens the scientific conclusion from the prior no-xp3 experiment.  Once lower moments are supplied by SAM, the present predictor needs a real turbulent-transport source; adding it restores amplitude, especially for r_t'^3.  But matching a peak is not the same as matching the vertical structure.  The poor same-time profile correlation means frozen two-Gaussian geometry alone does not determine the correct scalar-third-moment flux divergence."
    )
    report.heading("Recommendation")
    report.callout(
        "Do not tune c_K3 or nu3 to compensate for this yet.",
        "Keep the new budget terms permanently.  Next compare the diagnosed <w'r_t'^3> and <w'theta_l'^3> themselves against SAM before their divergence, and add the missing forcing covariance 3<x'^2 F_x'> only after microphysics/radiation can supply it consistently.  In parallel, evaluate whether a bounded PDF-geometry response based on the now-smooth LES-supplied lower moments improves flux placement.  The residual is no longer the uncertainty; the physical source and frozen-geometry assumptions are.",
        tone="warning",
    )
    report.heading("Reproducibility")
    report.paragraph(
        "The figures use immutable artifacts under agent_artifacts/xp3-budget-arm-20260724 and the checked-in ARM SAM archive, never output/.  The rendered PNGs, compact metrics JSON, and this exact plotting script are copied into the static bundle."
    )
    report.code_file(Path(__file__), language="python", caption="Reproduction script")
    (report.stage / "data").mkdir(exist_ok=True)
    shutil.copy2(payload_path, report.stage / "data" / payload_path.name)
    report.publish()


def main() -> int:
    runs, raw = _load()
    try:
        with tempfile.TemporaryDirectory(prefix="clubb_pdf10_xp3_budget_") as temp_dir:
            work = Path(temp_dir)
            metrics = _metrics(runs, raw)
            _plot_profiles(work, runs, raw)
            _plot_budget_profiles(work, runs["corrected"])
            _plot_time_height(work, runs["corrected"])
            _publish(work, metrics)
    finally:
        raw.close()
        for dataset in runs.values():
            dataset.close()
    print(f"Published doc/reports/{REPORT_ID}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
