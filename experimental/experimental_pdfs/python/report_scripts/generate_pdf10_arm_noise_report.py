#!/usr/bin/env python3
"""Publish a reproducible ARM PDF-10 scalar-third-moment noise audit.

This report intentionally reads the isolated SCM artifacts under
``agent_artifacts/``.  It never reads ``output/`` so later interactive runs
cannot silently change the scientific evidence in the published report.
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
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import xarray as xr

from dash_app.reports_tab.publisher import ReportBuilder


ARTIFACT_ROOT = ROOT / "agent_artifacts" / "xp3-noise-arm-20260724"
RAW_SAM = (
    ROOT
    / "input/les_and_clubb_benchmark_runs/sam_benchmark_runs/JULY_2017/ARM_96x96x110"
    / "GCSSARM_96x96x110_67m_40m_1s.nc"
)
REPORT_ID = "pdf10-arm-scalar-third-moment-noise-audit"
REPORT_TITLE = "PDF-10 ARM audit: where the scalar-third-moment noise enters"
SNAPSHOT_SECONDS = 69_780.0
HEIGHT_LIMIT_M = 4_380.0
# The ARM SAM archive records minutes from 11:30 UTC, whereas CLUBB SCM output
# uses seconds since midnight.  Preserve the physical clock when scoring.
SAM_CLOCK_OFFSET_SECONDS = 11.5 * 3600.0

RUNS = {
    "baseline": "baseline",
    "4x xp3 diffusion": "diffusion-4x",
    "30 s step": "half-timestep",
    "LES scalar 2nd moments": "les-xp2",
    "LES w moments": "les-wp23",
    "LES lower moments": "les-lower-moments",
}
COLORS = {
    "baseline": "#f97316",
    "4x xp3 diffusion": "#a78bfa",
    "30 s step": "#38bdf8",
    "LES scalar 2nd moments": "#ec4899",
    "LES w moments": "#22c55e",
    "LES lower moments": "#facc15",
}


def _revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def _field(ds: xr.Dataset, name: str) -> np.ndarray:
    """Return the single-column time-height field as a float NumPy array."""
    value = ds[name]
    # CLUBB uses ``col`` whereas archived SAM benchmark profiles retain
    # singleton x/y dimensions.  These artifacts are deliberately one-column
    # comparisons, so reduce any non-time/non-vertical singleton coordinate.
    vertical = "zm" if "zm" in value.dims else "zt" if "zt" in value.dims else "z"
    for dimension in tuple(value.dims):
        if dimension not in {"time", vertical}:
            value = value.isel({dimension: 0})
    return np.asarray(value.values, dtype=float)


def _coordinate(ds: xr.Dataset, name: str) -> np.ndarray:
    return np.asarray(ds[name].values, dtype=float)


def _z_for(ds: xr.Dataset, name: str) -> np.ndarray:
    return _coordinate(ds, "zm" if "zm" in ds[name].dims else "zt")


def _vertical_roughness(values: np.ndarray) -> float:
    """RMS three-point high-pass residual; smaller means vertically smoother."""
    if values.shape[1] < 3:
        return float("nan")
    # x_k - (x_{k-1} + 2 x_k + x_{k+1}) / 4: the departure from the
    # local three-point mean.  Retain the 0.5 center coefficient explicitly;
    # omitting it would measure curvature plus absolute amplitude.
    residue = 0.5 * values[:, 1:-1] - 0.25 * values[:, :-2] - 0.25 * values[:, 2:]
    return float(np.sqrt(np.nanmean(residue * residue)))


def _interp_raw(raw_values: np.ndarray, raw_time: np.ndarray, raw_z: np.ndarray,
                clubb_time: np.ndarray, clubb_z: np.ndarray) -> np.ndarray:
    """Interpolate raw SAM profile moments to the CLUBB time/height grid."""
    temporally = np.empty((len(clubb_time), len(raw_z)), dtype=float)
    for k in range(len(raw_z)):
        temporally[:, k] = np.interp(clubb_time, raw_time, raw_values[:, k])
    output = np.empty((len(clubb_time), len(clubb_z)), dtype=float)
    for t in range(len(clubb_time)):
        output[t, :] = np.interp(clubb_z, raw_z, temporally[t, :], left=np.nan, right=np.nan)
    return output


def _relative_rmse(values: np.ndarray, reference: np.ndarray) -> float:
    mask = np.isfinite(values) & np.isfinite(reference)
    numerator = np.sqrt(np.mean((values[mask] - reference[mask]) ** 2))
    denominator = np.sqrt(np.mean(reference[mask] ** 2))
    return float(numerator / denominator) if denominator > 0 else float("nan")


def _style_axis(axis: plt.Axes) -> None:
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
    figure.savefig(path, dpi=150, facecolor=figure.get_facecolor())
    plt.close(figure)


def _load_runs() -> dict[str, xr.Dataset]:
    data: dict[str, xr.Dataset] = {}
    for label, directory in RUNS.items():
        path = ARTIFACT_ROOT / directory / "arm_stats.nc"
        if not path.is_file():
            raise FileNotFoundError(f"required isolated SCM artifact is missing: {path}")
        data[label] = xr.open_dataset(path, decode_times=False)
    return data


def _metrics(runs: dict[str, xr.Dataset], raw: xr.Dataset) -> dict:
    baseline = runs["baseline"]
    raw_time = _coordinate(raw, "time") * 60.0 + SAM_CLOCK_OFFSET_SECONDS
    raw_z = _coordinate(raw, "z")
    result: dict[str, dict] = {}
    for label, ds in runs.items():
        time = _coordinate(ds, "time")
        row: dict[str, float] = {}
        for clubb_name, raw_name in (("rtp3", "RTP3"), ("thlp3", "THLP3"),
                                     ("rtp2", "RTP2"), ("thlp2", "THLP2"),
                                     ("wprtp", "WPRTP"), ("wpthlp", "WPTHLP"),
                                     ("wp3", "WP3")):
            values = _field(ds, clubb_name)
            z = _z_for(ds, clubb_name)
            valid_z = z <= HEIGHT_LIMIT_M
            reference = _interp_raw(_field(raw, raw_name), raw_time, raw_z, time, z)
            row[f"{clubb_name}_roughness"] = _vertical_roughness(values[:, valid_z])
            row[f"{clubb_name}_relative_rmse"] = _relative_rmse(
                values[:, valid_z], reference[:, valid_z]
            )
        result[label] = row

    # The SAM roughness is evaluated after interpolation to the baseline grid.
    sam_metrics: dict[str, float] = {}
    time = _coordinate(baseline, "time")
    for clubb_name, raw_name in (("rtp3", "RTP3"), ("thlp3", "THLP3"), ("wp3", "WP3")):
        z = _z_for(baseline, clubb_name)
        valid_z = z <= HEIGHT_LIMIT_M
        reference = _interp_raw(_field(raw, raw_name), raw_time, raw_z, time, z)
        sam_metrics[f"{clubb_name}_roughness"] = _vertical_roughness(reference[:, valid_z])
    result["raw SAM (interpolated)"] = sam_metrics
    return result


def _plot_sensitivity(work: Path, metrics: dict) -> None:
    labels = list(RUNS)
    short = ["base", "4× diff.", "30 s", "LES xp2", "LES wp23", "LES lower"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.8), facecolor="#0d1525")
    for axis, name, raw_name, title in (
        (axes[0], "rtp3", "RTP3", "r_t third-moment vertical roughness"),
        (axes[1], "thlp3", "THLP3", "theta_l third-moment vertical roughness"),
    ):
        raw_value = metrics["raw SAM (interpolated)"][f"{name}_roughness"]
        values = [metrics[label][f"{name}_roughness"] / raw_value for label in labels]
        axis.bar(np.arange(len(labels)), values, color=[COLORS[label] for label in labels])
        axis.axhline(1.0, color="#facc15", linewidth=1.8, linestyle="--", label="raw SAM")
        axis.set(xticks=np.arange(len(labels)), xticklabels=short, ylabel="roughness / SAM", title=title)
        axis.set_yscale("log")
        axis.legend(frameon=False, labelcolor="#e2e8f0", fontsize=8)
        _style_axis(axis)
    fig.suptitle("The third-moment noise follows the w-moment pathway, not scalar variance alone", color="#f8fafc", fontsize=14)
    _save(fig, work / "sensitivity-roughness.png")


def _plot_profiles(work: Path, runs: dict[str, xr.Dataset], raw: xr.Dataset) -> None:
    baseline = runs["baseline"]
    time = _coordinate(baseline, "time")
    ti = int(np.argmin(np.abs(time - SNAPSHOT_SECONDS)))
    raw_ti = int(np.argmin(np.abs(
        _coordinate(raw, "time") * 60.0 + SAM_CLOCK_OFFSET_SECONDS - time[ti]
    )))
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 6.6), facecolor="#0d1525", sharey=True)
    entries = (("wp3", "WP3", r"w'^3$ [m$^3$ s$^{-3}$]"),
               ("rtp3", "RTP3", r"r_t'^3$ [(kg kg$^{-1}$)$^3$]"),
               ("thlp3", "THLP3", r"theta_l'^3$ [K$^3$]"))
    for axis, (name, raw_name, xlabel) in zip(axes, entries):
        z = _z_for(baseline, name); mask = z <= 2800.0
        axis.plot(_field(raw, raw_name)[raw_ti], _coordinate(raw, "z"), color="#facc15", linewidth=2.0, label="raw SAM")
        for label in ("baseline", "LES scalar 2nd moments", "LES w moments", "LES lower moments"):
            ds = runs[label]
            axis.plot(_field(ds, name)[ti, mask], z[mask], color=COLORS[label], linewidth=1.55, label=label)
        axis.axvline(0.0, color="#64748b", linewidth=.8)
        axis.set(xlabel=xlabel, ylim=(0, 2800), title=name)
        _style_axis(axis)
    axes[0].set_ylabel("height [m]")
    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", ncol=3, frameon=False, labelcolor="#e2e8f0", fontsize=8)
    fig.suptitle(f"ARM snapshot at {time[ti] / 3600:.2f} h: forcing w moments removes the jagged component", color="#f8fafc", fontsize=14)
    fig.subplots_adjust(bottom=.16, top=.91, wspace=.24)
    fig.savefig(work / "moment-pathway-profiles.png", dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def _plot_skewness(work: Path, baseline: xr.Dataset) -> dict:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.5), facecolor="#0d1525")
    output: dict[str, float] = {}
    for axis, variance, skewness, threshold, title in (
        (axes[0], "rtp2", "Skrt_zt", 1.0e-8, r"r_t: |Skrt| grows where r_t variance is floored"),
        (axes[1], "thlp2", "Skthl_zt", 1.0e-2, r"theta_l: |Skthl| grows where theta_l variance is floored"),
    ):
        # rtp2/thlp2 are zm; simple centered interpolation gives a diagnostic zt value.
        var_m = _field(baseline, variance)
        var_t = 0.5 * (var_m[:, :-1] + var_m[:, 1:])
        sk = np.abs(_field(baseline, skewness))
        mask = np.isfinite(var_t) & np.isfinite(sk) & (var_t > 0) & (sk > 0)
        axis.hexbin(np.log10(var_t[mask]), np.log10(sk[mask]), gridsize=48, mincnt=1, cmap="magma")
        axis.axvline(np.log10(threshold), color="#facc15", linestyle="--", linewidth=1.4)
        axis.set(xlabel=f"log10({variance})", ylabel=f"log10(|{skewness}|)", title=title)
        _style_axis(axis)
        output[f"fraction_{variance}_below_threshold"] = float(np.mean(var_t < threshold))
        visible = sk[var_t >= threshold]
        output[f"{skewness}_rms_above_threshold"] = float(np.sqrt(np.mean(visible * visible)))
    fig.suptitle("Huge normalized-skewness values are often a denominator diagnostic, not a usable resolved signal", color="#f8fafc", fontsize=14)
    _save(fig, work / "skewness-denominator.png")
    return output


def _plot_wp3_budgets(work: Path, baseline: xr.Dataset) -> dict:
    time = _coordinate(baseline, "time")
    ti = int(np.argmin(np.abs(time - SNAPSHOT_SECONDS)))
    zt = _coordinate(baseline, "zt")
    zm = _coordinate(baseline, "zm")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 6.2), facecolor="#0d1525")
    for name, color in (("wp3_ta", "#f97316"), ("wp3_dp1", "#38bdf8"),
                        ("wp3_bp1", "#a78bfa"), ("wp3_cl", "#22c55e")):
        axes[0].plot(_field(baseline, name)[ti], zt, color=color, linewidth=1.6, label=name)
    axes[0].axvline(0.0, color="#64748b", linewidth=.8)
    axes[0].set(xlabel=r"tendency [m$^3$ s$^{-4}$]", ylabel="height [m]", ylim=(0, 2800), title="wp3 budget terms")
    axes[0].legend(frameon=False, labelcolor="#e2e8f0", fontsize=8)
    _style_axis(axes[0])
    coefficient = _field(baseline, "coef_wp4_implicit")
    axes[1].plot(coefficient[ti], zm, color="#facc15", linewidth=1.9, label="PDF-10 frozen wp4/wp2²")
    axes[1].set(xlabel="coef_wp4_implicit [-]", ylabel="height [m]", ylim=(0, 2800), title="PDF fourth-moment response")
    axes[1].legend(frameon=False, labelcolor="#e2e8f0", fontsize=8)
    _style_axis(axes[1])
    fig.suptitle(f"wp3 transport and dissipation are the next terms to close at {time[ti] / 3600:.2f} h", color="#f8fafc", fontsize=14)
    _save(fig, work / "wp3-budget-and-response.png")
    return {
        "snapshot_seconds": float(time[ti]),
        "coefficient_max": float(np.nanmax(coefficient)),
        "coefficient_p99": float(np.nanpercentile(coefficient, 99.0)),
        "wp3_ta_roughness": _vertical_roughness(_field(baseline, "wp3_ta")),
        "wp3_dp1_roughness": _vertical_roughness(_field(baseline, "wp3_dp1")),
    }


def _publish(work: Path, metrics: dict, skewness: dict, budget: dict) -> None:
    baseline = metrics["baseline"]
    raw = metrics["raw SAM (interpolated)"]
    wp23 = metrics["LES w moments"]
    diffusion = metrics["4x xp3 diffusion"]
    payload = {"runs": RUNS, "metrics": metrics, "skewness": skewness, "wp3_budget": budget}
    data_path = work / "noise-audit-data.json"
    data_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    report = ReportBuilder(
        REPORT_ID,
        REPORT_TITLE,
        summary="A controlled ARM SCM audit finds that PDF-10 r_t'^3 and theta_l'^3 noise is injected chiefly through the free-running w'^2/w'^3/PDF-response pathway.  Stronger scalar-third-moment diffusion can hide it, but is not a clean fix; variance-normalized skewness also becomes uninformative over the widespread scalar-variance floor.",
        tags=("PDF-10", "ARM", "scalar third moments", "noise audit", "LES advance"),
        replace=True,
        source_revision=_revision(),
    )
    report.heading("Bottom line")
    report.callout(
        "The leading source is upstream of the new r_t'^3/theta_l'^3 solver.",
        "Forcing only the SAM w moments (wp2 and wp3) reduces the raw r_t third-moment vertical roughness by about %.0fx and the theta_l third-moment roughness by about %.0fx.  Forcing scalar second moments alone does not.  The immediate investigation target is therefore the free-running wp3 budget and the PDF fourth-moment response that enters its semi-implicit turbulent-advection operator." % (
            baseline["rtp3_roughness"] / wp23["rtp3_roughness"],
            baseline["thlp3_roughness"] / wp23["thlp3_roughness"],
        ),
        tone="warning",
    )
    report.metrics([
        ("r_t'^3 roughness", f"{baseline['rtp3_roughness'] / raw['rtp3_roughness']:.1f}× SAM", "baseline vertical high-pass RMS"),
        ("theta_l'^3 roughness", f"{baseline['thlp3_roughness'] / raw['thlp3_roughness']:.1f}× SAM", "baseline vertical high-pass RMS"),
        ("w-moment forcing", f"{baseline['rtp3_roughness'] / wp23['rtp3_roughness']:.1f}× smoother", "r_t'^3 after forcing wp2/wp3 only"),
        ("PDF wp4 response", f"p99 {budget['coefficient_p99']:.1f}", f"maximum {budget['coefficient_max']:.1f}; inspect, do not cap blindly"),
    ])
    report.heading("Experiment design")
    report.paragraph(
        "Each comparison is an isolated ARM SCM run in agent_artifacts/xp3-noise-arm-20260724, not the mutable output directory.  All runs use PDF-10 and the same current PDF-10 settings; they differ by one numerical sensitivity or one LES-advance group.  The raw SAM source is interpolated only for common time/height scoring.  Vertical roughness is the RMS of a three-point high-pass residual, so it measures grid-scale jaggedness rather than physical profile amplitude."
    )
    report.code(
        "baseline: same PDF-10 parameters as the current ARM configuration\n"
        "4x xp3 diffusion: c_K3=0.10, nu3=5 (baseline: 0.025, 1)\n"
        "30 s step: dt_main=30 s (baseline: 60 s)\n"
        "LES scalar 2nd moments: -les_advance xp2\n"
        "LES w moments: -les_advance wp23\n"
        "LES lower moments: -les_advance wp23,xm,xp2,windm",
        language="text", caption="Controlled sensitivity matrix",
    )
    report.heading("Evidence: the w-moment pathway controls the raw third-moment roughness")
    report.figure(work / "sensitivity-roughness.png", caption="Vertical roughness relative to raw SAM.  Fourfold scalar-third-moment diffusion suppresses roughness, but forcing wp2/wp3 removes it much more selectively.  Scalar second-moment forcing alone barely changes the raw scalar-third-moment roughness.")
    report.figure(work / "moment-pathway-profiles.png", caption="At an ARM 19.38 h snapshot, the wp23-forced run follows the smooth raw-SAM wp3 structure and the r_t/theta_l third moments lose their jagged profile.  Exact scalar second moments alone do not accomplish this.  The remaining amplitude bias is separate from the noise source.")
    report.heading("What this rules out")
    report.paragraph(
        "Halving the main timestep reduces r_t'^3 roughness modestly but leaves theta_l'^3 essentially unchanged; this is not a simple 60 s time-step instability.  The new semi-implicit scalar-third-moment solve also becomes smooth when wp2/wp3 are supplied from SAM, so it is not by itself generating the grid-scale noise.  Finally, using only LES scalar second moments makes rtp2 and thlp2 smooth but leaves rtp3/thlp3 rough, ruling out scalar-variance noise as the dominant direct injector."
    )
    report.heading("Why more c_K3 / nu3 is not the recommended fix")
    report.paragraph(
        "The 4x diffusion sensitivity lowers r_t'^3 and theta_l'^3 roughness to %.2fx and %.2fx of SAM, respectively.  That is useful evidence that the high-frequency mode is diffusive, but it also substantially changes the evolved scalar-third-moment field and can create enormous normalized skewness where r_t variance is nearly zero.  It should be retained as a guardrail sensitivity, not promoted to a tuned default until the upstream wp3/PDF response is understood." % (
            diffusion["rtp3_roughness"] / raw["rtp3_roughness"], diffusion["thlp3_roughness"] / raw["thlp3_roughness"],
        )
    )
    report.heading("A second, separate issue: normalized skewness is often undefined in practice")
    report.figure(work / "skewness-denominator.png", caption="The enormous Skrt_zt and Skthl_zt values cluster where rtp2/thlp2 is at or near its numerical floor.  These plots do not prove that the raw third moments are correct; they show that a ratio by a tiny variance cannot be used as a reliable noise score there.")
    report.paragraph(
        "In this baseline, %.1f%% of r_t variance samples fall below 1e-8 (kg/kg)^2 and %.1f%% of theta_l variance samples fall below 1e-2 K^2.  Those are analysis-screen thresholds, not proposed closure constants.  Before changing Skx_module or its denominator tolerance, audit every downstream use.  The safe immediate change is in diagnostics/scoring: display raw xp3 and label or mask normalized skewness below a stated variance threshold." % (
            100.0 * skewness["fraction_rtp2_below_threshold"], 100.0 * skewness["fraction_thlp2_below_threshold"],
        )
    )
    report.heading("The concrete upstream hypothesis: wp3 transport, dissipation, and the frozen PDF response")
    report.figure(work / "wp3-budget-and-response.png", caption="The free-running wp3 budget has comparatively rough turbulent-advection and dissipation terms.  The semi-implicit operator uses the diagnosed PDF ratio wp4/wp2²; it is positive as required but has a high upper tail.  This is a plausible gain pathway, not proof that the coefficient itself is erroneous.")
    report.paragraph(
        "PDF-10 diagnoses coef_wp4_implicit from its two-Gaussian geometry, and advance_wp2_wp3_module.F90 uses it in the wp3 turbulent-advection band operator.  The audit therefore points to the coupled wp3/PDF geometry loop.  A large coefficient can be physically meaningful for intermittent vertical velocity, but any level-to-level jump in it can also alter the implicit operator strongly.  Its vertical structure, its relation to wp2 and mixture fraction, and exact wp3 budget closure need to be scored against raw SAM before adding a cap or smoothing it."
    )
    report.heading("Recommendations, in order")
    report.code(
        "1. Instrument the actual PDF-10 xp3 budget in standard_stats.in.\n"
        "   Register rtp3_ta/thlp3_ta plus production, accumulation, mean-advection,\n"
        "   diffusion, and damping contributions emitted by advance_xp3_pdf10.\n\n"
        "2. Audit wp3 before tuning scalar-third-moment diffusion.\n"
        "   Score wp3_ta, wp3_dp1, wp3_bp1, wp3_cl and coef_wp4_implicit versus\n"
        "   SAM wp3 at the same output cadence; test whether a frozen/smoothed\n"
        "   response changes the noisy high-pass component without moving the mean.\n\n"
        "3. Keep the 30 s and 4x-diffusion cases as regression sensitivities only.\n"
        "   The former diagnoses time-step sensitivity; the latter gives an upper\n"
        "   bound on what pure scalar-xp3 damping can hide.\n\n"
        "4. Make plot/tuner diagnostics variance-aware now.\n"
        "   Score Skrt/Skthl only above explicit variance thresholds and always show\n"
        "   raw rtp3/thlp3 beside them.  Do not change closure denominator constants\n"
        "   until consumers of Skx are audited.\n\n"
        "5. After wp3 is smoother, reassess missing xp3 physics.\n"
        "   The wp23-forced solution is smooth but still under-amplitude, indicating\n"
        "   the scalar-third-moment equation needs better sources/transport after the\n"
        "   upstream injector is fixed—not simply more diffusion.",
        language="text", caption="Recommended next development sequence",
    )
    report.heading("Known budget-output gap")
    report.paragraph(
        "The current PDF-10 scalar-third-moment routine emits its new tendency labels (for example rtp3_ta and thlp3_ta), but standard_stats.in does not yet register them.  The legacy rtp3_bt/thlp3_bt outputs are therefore not a valid closure check for this path.  Adding the explicit term labels is a small observability change and should precede scientific retuning."
    )
    report.heading("Reproducibility")
    report.paragraph(
        "All figures derive from the retained, dedicated SCM outputs listed in the embedded JSON.  The report script is included verbatim.  It does not use output/arm_stats.nc, so rerunning ARM later cannot revise this evidence.  The static report bundle contains its rendered PNGs and needs no live calculation to view."
    )
    report.code_file(Path(__file__), language="python", caption="Reproduction script")
    (report.stage / "data").mkdir(exist_ok=True)
    shutil.copy2(data_path, report.stage / "data" / data_path.name)
    report.publish()


def main() -> int:
    if not RAW_SAM.is_file():
        raise FileNotFoundError(f"raw ARM SAM source is unavailable: {RAW_SAM}")
    runs = _load_runs()
    raw = xr.open_dataset(RAW_SAM, decode_times=False)
    try:
        with tempfile.TemporaryDirectory(prefix="clubb_pdf10_noise_report_") as directory:
            work = Path(directory)
            metrics = _metrics(runs, raw)
            _plot_sensitivity(work, metrics)
            _plot_profiles(work, runs, raw)
            skewness = _plot_skewness(work, runs["baseline"])
            budget = _plot_wp3_budgets(work, runs["baseline"])
            _publish(work, metrics, skewness, budget)
    finally:
        raw.close()
        for ds in runs.values():
            ds.close()
    print(f"Published doc/reports/{REPORT_ID}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
