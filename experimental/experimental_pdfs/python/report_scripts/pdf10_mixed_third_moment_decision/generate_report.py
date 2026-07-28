#!/usr/bin/env python3
"""Evaluate mixed third moments as independent PDF-10 center information.

This script compares the current full-LES-advance ARM PDF-10 diagnostic
moments with the matching raw SAM profile moments.  It publishes a static
decision report; viewing the report never reruns this code.
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


REPORT_ID = "pdf10-mixed-third-moment-decision"
REPORT_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPORT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_pdf10_mixed_moment_report_matplotlib")

from dash_app.reports_tab.publisher import ReportBuilder


RAW_PATH = REPO_ROOT / "input/les_and_clubb_benchmark_runs/sam_benchmark_runs/JULY_2017/ARM_96x96x110/GCSSARM_96x96x110_67m_40m_1s.nc"
STATS_PATH = REPO_ROOT / "output/arm_stats.nc"
EXAMPLE_MODEL_SECONDS = 71_460.0  # 19.85 h; raw SAM minute 501.
EXAMPLE_HEIGHT_M = 540.0
MAX_ANALYSIS_HEIGHT_M = 2_600.0


TERMS = (
    {
        "raw": "WPRTP2", "clubb": "wprtp2", "label": r"$\overline{w'r_t'^2}$",
        "description": "vertical-motion sign of large-magnitude moisture anomalies",
    },
    {
        "raw": "WP2RTP", "clubb": "wp2rtp", "label": r"$\overline{w'^2r_t'}$",
        "description": "moisture sign carried by energetic vertical motion",
    },
    {
        "raw": "WPRTPTHLP", "clubb": "wprtpthlp", "label": r"$\overline{w'r_t'\theta_l'}$",
        "description": "vertical-motion partitioning of the thermodynamic covariance",
    },
)


def _revision() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, text=True).strip() + " (working tree)"
    except (OSError, subprocess.CalledProcessError):
        return "working tree"


def _profile(dataset, name):
    values = np.asarray(dataset[name].values, dtype=float)
    return values[:, :, 0, 0] if values.ndim == 4 else values[:, :, 0]


def _normalized_terms(raw, clubb):
    """Return raw/PDF normalized third moments on the common 40 m z grid."""
    wp2 = _profile(clubb, "wp2")[:, :raw.sizes["z"]]
    rtp2 = _profile(clubb, "rtp2")[:, :raw.sizes["z"]]
    thlp2 = _profile(clubb, "thlp2")[:, :raw.sizes["z"]]
    scales = (
        np.sqrt(np.maximum(wp2, 0.0)) * rtp2,
        wp2 * np.sqrt(np.maximum(rtp2, 0.0)),
        np.sqrt(np.maximum(wp2 * rtp2 * thlp2, 0.0)),
    )
    rows = []
    for meta, scale in zip(TERMS, scales):
        source = _profile(raw, meta["raw"])
        modeled = _profile(clubb, meta["clubb"])[:, :raw.sizes["z"]]
        rows.append({**meta, "raw_values": source, "pdf_values": modeled, "scale": scale})
    return rows


def _score(row, z):
    raw = row["raw_values"]
    pdf = row["pdf_values"]
    scale = row["scale"]
    admissible = (
        (z[None, :] <= MAX_ANALYSIS_HEIGHT_M)
        & np.isfinite(raw)
        & np.isfinite(pdf)
        & np.isfinite(scale)
        & (scale > 1.0e-20)
    )
    scale_floor = np.nanpercentile(scale[admissible], 20.0)
    mask = admissible & (scale > scale_floor)
    target = raw[mask] / scale[mask]
    prediction = pdf[mask] / scale[mask]
    strong = np.abs(target) > 0.08
    return {
        "states": int(mask.sum()),
        "correlation": float(np.corrcoef(target, prediction)[0, 1]),
        "normalized_rmse": float(np.sqrt(np.mean((prediction - target) ** 2))),
        "normalized_bias": float(np.mean(prediction - target)),
        "strong_signal_sign_agreement": float(np.mean(np.sign(target[strong]) == np.sign(prediction[strong]))),
        "strong_signal_states": int(strong.sum()),
    }


def _render_scatter(work_dir, rows, z):
    from matplotlib import pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(16.8, 5.5), dpi=170, facecolor="#0d1525")
    scores = []
    for axis, row in zip(axes, rows):
        score = _score(row, z)
        scores.append({"name": row["raw"], **score})
        raw, pdf, scale = row["raw_values"], row["pdf_values"], row["scale"]
        mask = (z[None, :] <= MAX_ANALYSIS_HEIGHT_M) & np.isfinite(raw) & np.isfinite(pdf) & (scale > 1.0e-20)
        floor = np.nanpercentile(scale[mask], 20.0)
        mask &= scale > floor
        target, prediction = raw[mask] / scale[mask], pdf[mask] / scale[mask]
        height = np.broadcast_to(z[None, :], raw.shape)[mask]
        limit = max(1.0, float(np.percentile(np.abs(np.concatenate((target, prediction))), 99.3)))
        dots = axis.scatter(target, prediction, c=height, s=3.2, cmap="viridis", alpha=0.45, rasterized=True)
        axis.plot((-limit, limit), (-limit, limit), color="#facc15", linewidth=1.4, linestyle="--")
        axis.axhline(0.0, color="#64748b", linewidth=.7); axis.axvline(0.0, color="#64748b", linewidth=.7)
        axis.set_xlim(-limit, limit); axis.set_ylim(-limit, limit); axis.set_aspect("equal", adjustable="box")
        axis.set_facecolor("#111827")
        axis.set_title(f"{row['label']}\nr={score['correlation']:.2f}, normalized RMSE={score['normalized_rmse']:.2f}", color="#edf2f7", fontsize=10)
        axis.set_xlabel("raw SAM normalized moment", color="#dbe5f1", fontsize=9)
        axis.set_ylabel("PDF-10 diagnosed normalized moment", color="#dbe5f1", fontsize=9)
        axis.tick_params(colors="#cbd5e1", labelsize=8)
        for spine in axis.spines.values(): spine.set_color("#40516a")
        colorbar = figure.colorbar(dots, ax=axis, pad=.015, shrink=.83)
        colorbar.set_label("height [m]", color="#cbd5e1", fontsize=8); colorbar.ax.tick_params(colors="#cbd5e1", labelsize=7)
    figure.suptitle("Mixed third moments: raw SAM targets versus PDF-10 diagnoses", color="#f8fafc", fontsize=15, y=.99)
    figure.text(.5, .015, "All panels use the current ARM full-LES-advance run.  Inputs through pure scalar skewness are forced; these mixed moments remain a consequence of the PDF geometry.", ha="center", color="#b6c5d8", fontsize=9)
    figure.tight_layout(rect=(0, .05, 1, .93))
    path = work_dir / "mixed-third-moment-scatter.png"; figure.savefig(path, facecolor=figure.get_facecolor()); plt.close(figure)
    return scores


def _render_profile_and_example(work_dir, rows, raw, clubb):
    from matplotlib import pyplot as plt

    z = np.asarray(raw.z.values, dtype=float)
    raw_index = int(np.argmin(np.abs(np.asarray(raw.time.values, dtype=float) - 501.0)))
    clubb_index = int(np.argmin(np.abs(np.asarray(clubb.time.values, dtype=float) - EXAMPLE_MODEL_SECONDS)))
    height_index = int(np.argmin(np.abs(z - EXAMPLE_HEIGHT_M)))
    figure, axes = plt.subplots(1, 2, figsize=(15.2, 6.2), dpi=170, facecolor="#0d1525")
    colors = ("#38bdf8", "#f97316", "#c084fc")
    example = []
    for color, row in zip(colors, rows):
        source = row["raw_values"][raw_index] / row["scale"][clubb_index]
        modeled = row["pdf_values"][clubb_index] / row["scale"][clubb_index]
        axes[0].plot(source, z, color=color, linewidth=2.2, label=f"SAM {row['label']}")
        axes[0].plot(modeled, z, color=color, linewidth=1.6, linestyle="--", label=f"PDF-10 {row['label']}")
        example.append({
            "name": row["raw"], "sam_normalized": float(source[height_index]),
            "pdf10_normalized": float(modeled[height_index]),
            "sam_raw": float(row["raw_values"][raw_index, height_index]),
            "pdf10_raw": float(row["pdf_values"][clubb_index, height_index]),
        })
    axes[0].axhline(EXAMPLE_HEIGHT_M, color="#facc15", linestyle="--", linewidth=1.2)
    axes[0].axvline(0.0, color="#64748b", linewidth=.8)
    axes[0].set(xlabel="normalized mixed third moment", ylabel="height [m]", ylim=(0, 2_600))
    axes[0].set_title("19.85 h profile: raw SAM solid, PDF-10 dashed", color="#edf2f7", fontsize=11)
    axes[0].legend(fontsize=8, loc="lower left", frameon=False, labelcolor="#e2e8f0")
    labels = [row["label"] for row in rows]
    raw_values = [item["sam_normalized"] for item in example]
    pdf_values = [item["pdf10_normalized"] for item in example]
    x = np.arange(len(labels)); width=.36
    axes[1].bar(x-width/2, raw_values, width, color="#facc15", label="raw SAM")
    axes[1].bar(x+width/2, pdf_values, width, color="#ef4444", label="PDF-10")
    axes[1].axhline(0.0, color="#94a3b8", linewidth=.8)
    axes[1].set_xticks(x, labels); axes[1].set_ylabel("normalized moment")
    axes[1].set_title("The low-level V example: 501 min / 540 m", color="#edf2f7", fontsize=11)
    axes[1].legend(frameon=False, labelcolor="#e2e8f0")
    for axis in axes:
        axis.set_facecolor("#111827"); axis.tick_params(colors="#cbd5e1", labelsize=9)
        axis.xaxis.label.set_color("#dbe5f1"); axis.yaxis.label.set_color("#dbe5f1")
        for spine in axis.spines.values(): spine.set_color("#40516a")
    figure.suptitle("Mixed moments see the V-type tail partition that marginal skewness loses", color="#f8fafc", fontsize=15, y=.99)
    figure.tight_layout(rect=(0, 0, 1, .94))
    path = work_dir / "mixed-third-moment-v-example.png"; figure.savefig(path, facecolor=figure.get_facecolor()); plt.close(figure)
    return example


def _publish(work_dir, scores, example):
    payload = {
        "raw_sam_profile": str(RAW_PATH), "clubb_output": str(STATS_PATH),
        "analysis_height_cap_m": MAX_ANALYSIS_HEIGHT_M,
        "example": {"model_time_seconds": EXAMPLE_MODEL_SECONDS, "raw_minute": 501.0, "height_m": EXAMPLE_HEIGHT_M},
        "scores": scores, "v_example": example,
    }
    data_path = work_dir / "mixed-third-moment-evaluation.json"
    data_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    score = {row["name"]: row for row in scores}
    report = ReportBuilder(
        REPORT_ID,
        "PDF-10 center placement: do mixed third moments justify a new prognostic closure?",
        summary="ARM full-LES-advance evidence says mixed third moments contain the missing conditional-tail information, but does not yet justify prognosing them: first isolate their payoff with a non-circular predictor experiment.",
        tags=("PDF-10", "ARM", "SAM", "mixed third moments", "prognostic moments", "LES advance"),
        replace=True, source_revision=_revision(),
    )
    report.heading("Decision")
    report.callout(
        "Do not add a new prognostic moment yet.",
        "The mixed third moments clearly expose information that PDF-10 does not use, and the current PDF diagnoses their V-case magnitudes poorly.  But one moment does not uniquely specify the V, while a genuine prognostic third-moment equation needs unclosed fourth-order transport and new stability work.  First perform a controlled, LES-forced predictor experiment using the existing SAM mixed-moment truth; only then decide whether its benefit warrants a new prognostic equation.",
        tone="warning",
    )
    report.metrics([
        ("Best pattern correlation", f"{score['WP2RTP']['correlation']:.2f}", "raw versus PDF-10 normalized WP2RTP below 2.6 km"),
        ("V: WP2RTP", f"{example[1]['sam_normalized']:+.2f} → {example[1]['pdf10_normalized']:+.2f}", "raw SAM versus PDF-10 at 501 min / 540 m"),
        ("V: WPRTPTHLP", f"{example[2]['sam_normalized']:+.2f} → {example[2]['pdf10_normalized']:+.2f}", "the trivariate thermodynamic partition is severely muted"),
        ("Recommendation", "predictor first", "not an exact additional-moment fit"),
    ])
    report.heading("What the existing terms add")
    report.paragraph(
        "Pure skewnesses say separately that w has a positive tail, r_t a negative tail, and theta_l a positive tail.  They do not say whether those tails are in the same parcels.  WP2RTP asks whether energetic vertical motion is moist or dry; WPRTP2 asks whether strong moisture departures rise or sink; WPRTPTHLP asks how vertical motion partitions the r_t–theta_l covariance.  Together they provide precisely the missing conditional information behind the V."
    )
    report.heading("Evidence: these are independent SAM targets, not merely pretty diagnostics")
    report.paragraph(
        "The raw SAM profile file contains all three quantities.  PDF-10 currently diagnoses them after selecting its components, so their disagreement is a closure error rather than an input mismatch.  The scatter comparison uses variance-normalized moments and excludes small-variance upper levels.  Pattern correlations are meaningful, but the systematic amplitude errors show that the diagnosed values cannot be used as an independent predictor without circularity."
    )
    report.figure(work_dir / "mixed-third-moment-scatter.png", caption="Raw SAM versus PDF-10-diagnosed normalized mixed third moments below 2.6 km.  The 1:1 line is gold.  PDF-10 retains some broad pattern but has substantial normalized error and bias, especially for the tail partition moments needed by the V.")
    report.figure(work_dir / "mixed-third-moment-v-example.png", caption="At the low-level V, all three SAM terms identify organized thermodynamic/vertical-motion structure.  PDF-10 greatly underdiagnoses WP2RTP and WPRTPTHLP while overdiagnosing WPRTP2, consistent with its misplaced/broad component allocation.")
    report.heading("Recommended next experiment: no prognostic code")
    report.paragraph(
        "Extend the LES-advance prototype only enough to read and temporarily supply SAM WP2RTP plus WPRTPTHLP to PDF-10 before center placement.  Use them as bounded compatibility predictors, not exact constraints: for example, shrink or rotate a marginal-skewness center direction when it conflicts with the sign of energetic-motion moisture and thermodynamic partitioning.  Re-run ARM, BOMEX, and at least one stratocumulus case, scoring component-center placement, wprcp, and the mixed moments themselves.  This directly measures the value of the missing information without claiming that free-running CLUBB can yet predict it."
    )
    report.heading("If that experiment succeeds: where a prognostic implementation belongs")
    report.paragraph(
        "Start with WP2RTP, optionally paired with WPRTPTHLP; do not add all five mixed moments at once.  WP2RTP already has a core-state/API route and is consumed by the mean/flux turbulent-advection path, but PDF closure currently overwrites it diagnostically.  A true independent prediction would require a new third-moment tendency owner, then an explicit decision about when PDF-10 may consume it rather than overwrite it.  WPRTP2 and WPRTPTHLP are presently PDF outputs and would require additional state plumbing as well."
    )
    report.code(
        """Implementation map if (and only if) the LES-forced experiment demonstrates robust benefit\n\n1. State + API\n   - retain/add independent mixed-moment state in clubb_driver.F90 and clubb_api_module.F90\n   - preserve it through initialization, restart, and grid_adaptation_module.F90\n\n2. Tendency owner\n   - add a dedicated advance_wp2xp / mixed-third-moment module\n   - include production, pressure/redistribution, dissipation, and source terms\n   - close turbulent transport of the third moment (a fourth-order flux); this is the hard part\n\n3. PDF interaction\n   - in pdf_closure_module.F90, do not overwrite the independently advanced input\n   - add a bounded PDF-10 compatibility/predictor use before center allocation\n   - continue diagnosing other PDF outputs from the final components\n\n4. LES-advance validation plumbing\n   - add normalized fields and source mappings in utilities/benchmark_converter.py\n   - add an explicit LES-advance group/flags in les_advancer_module.F90\n   - wire before/after endpoint forcing in advance_clubb_core_module.F90\n\n5. Budgets and safety\n   - add tendency-budget stats, realizability limits, and restart tests\n   - validate free-running stability before treating LES-advance agreement as closure success""",
        language="text", caption="A new prognostic mixed third moment is materially more invasive than adding a PDF-10 parameter.")
    report.heading("Why explicit second-moment turbulent advection is not the answer")
    report.paragraph(
        "The existing l_explicit_turbulent_adv_xpyp pathway governs turbulent advection of scalar second moments.  It does not provide a prognostic equation for a mixed third moment.  Advancing WP2RTP or WPRTPTHLP would still require their own tendencies and a closure for their fourth-order turbulent transport.  That is why the predictor experiment is the appropriate gate."
    )
    report.heading("Reproducibility")
    report.paragraph("Source profiles are the ARM SAM file GCSSARM_96x96x110_67m_40m_1s.nc; modeled diagnostics are output/arm_stats.nc from the current PDF-10 full-LES-advance run.  The report’s compact JSON records normalization, scores, and the low-level V values.")
    report.code_file(Path(__file__), language="python", caption="Reproduction script")
    (report.stage / "data").mkdir(exist_ok=True)
    shutil.copy2(data_path, report.stage / "data" / data_path.name)
    report.publish()


def main() -> int:
    if not RAW_PATH.is_file() or not STATS_PATH.is_file():
        raise FileNotFoundError("ARM raw SAM or current CLUBB stats output is unavailable.")
    with tempfile.TemporaryDirectory(prefix="clubb_mixed_third_report_") as directory:
        work_dir = Path(directory)
        raw = xr.open_dataset(RAW_PATH, decode_times=False)
        clubb = xr.open_dataset(STATS_PATH, decode_times=False)
        try:
            rows = _normalized_terms(raw, clubb)
            scores = _render_scatter(work_dir, rows, np.asarray(raw.z.values, dtype=float))
            example = _render_profile_and_example(work_dir, rows, raw, clubb)
            _publish(work_dir, scores, example)
        finally:
            raw.close(); clubb.close()
    print(f"Published doc/reports/{REPORT_ID}/report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
