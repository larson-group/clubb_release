#!/usr/bin/env python3
"""Generate the static directional-Lscale PDF-10 coherence report.

Run from repository root with:
  PYTHONPATH=.:.venv-dash/lib/python3.14/site-packages .venv-dash/bin/python \
    doc/reports/directional-lscale-pdf-coherence/snippets/generate_report.py
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dash_app.misc_tab.reports import trivariate_transport_sam_lab as base
from dash_app.misc_tab.sam_lscale import compute_sam_lscale_profile
from dash_app.misc_tab.transport_2g_prototype import diagnose_transport_2g_from_moments
from dash_app.misc_tab.vertical_coherence import (
    VerticalCoherenceSettings,
    apply_local_vertical_coherence,
    standardized_displacement,
)
from dash_app.reports_tab.publisher import ReportBuilder


REPORT_ID = "directional-lscale-pdf-coherence"
TIME_SECONDS = 52200
TARGET_HEIGHT_M = 540.0
FIGURE_DIR = Path("/tmp") / "clubb_directional_lscale_report"


def _ellipse(mean: np.ndarray, covariance: np.ndarray, radius: float = 2.0):
    scale = np.diag((1000.0, 1.0))
    cov = scale @ covariance[np.ix_((1, 2), (1, 2))] @ scale
    values, vectors = np.linalg.eigh(0.5 * (cov + cov.T))
    angles = np.linspace(0.0, 2.0 * np.pi, 181)
    points = np.array((mean[1] * 1000.0, mean[2]))[:, None] + radius * vectors @ np.diag(np.sqrt(np.maximum(values, 0.0))) @ np.vstack(
        (np.cos(angles), np.sin(angles))
    )
    return points


def _write_rt_thl_figure(snapshot, local, coherent, path: Path):
    counts, xedges, yedges = np.histogram2d(
        snapshot.samples[:, 1] * 1000.0, snapshot.samples[:, 2], bins=72
    )
    figure, axis = plt.subplots(figsize=(9.8, 6.4), constrained_layout=True)
    mesh = axis.pcolormesh(xedges, yedges, np.log1p(counts.T), cmap="cividis", shading="auto")
    figure.colorbar(mesh, ax=axis, label="log(1 + raw SAM count)")
    colors = ("#fbbf24", "#22d3ee")
    for result, style, label in ((local, "--", "local"), (coherent, "-", "directional Lscale prior")):
        for number, (mean, covariance, color) in enumerate(
            zip(result.mixture.means, result.mixture.covariances, colors), start=1
        ):
            curve = _ellipse(mean, covariance)
            axis.plot(
                curve[0], curve[1], style, color=color, linewidth=2.2,
                label=f"{label} G{number}" if number == 1 else None,
            )
            axis.scatter(mean[1] * 1000.0, mean[2], marker="x", color=color, s=55, zorder=5)
    axis.set(
        title=f"ARM raw SAM rₜ–θₗ at {snapshot.elapsed_minutes:.0f} min, {snapshot.height_m:.0f} m",
        xlabel="Total water, rₜ [g/kg]",
        ylabel="Liquid-water potential temperature, θₗ [K]",
    )
    axis.legend(frameon=False, loc="best")
    figure.savefig(path, dpi=180, facecolor="#0d1525")
    plt.close(figure)


def _write_lscale_figure(profile, target_height: float, path: Path):
    z = profile.height_m
    figure, axis = plt.subplots(figsize=(8.8, 7.0), constrained_layout=True)
    step = max(1, len(z) // 26)
    for values, color, sign, label in (
        (profile.lscale_up_m, "#38bdf8", 1.0, "Lscale_up"),
        (profile.lscale_down_m, "#fb7185", -1.0, "Lscale_down"),
    ):
        for index in range(0, len(z), step):
            start, arrival = z[index], z[index] + sign * max(values[index], 0.0)
            axis.plot((start, start), (start, arrival), color=color, alpha=0.70, linewidth=1.2)
            axis.scatter(start, start, s=12, facecolors="none", edgecolors=color)
            axis.scatter(start, arrival, s=35, marker="^" if sign > 0 else "v", color=color)
    axis.plot((z[0], z[-1]), (z[0], z[-1]), ":", color="#94a3b8", linewidth=1.1)
    axis.axhline(target_height, color="#fbbf24", linestyle="--", linewidth=1.3)
    axis.set(
        title="Native calc_Lscale_directly reach from one raw-SAM ARM column",
        xlabel="Parcel launch altitude [m]",
        ylabel="Estimated arrival altitude [m]",
    )
    axis.legend(
        handles=[
            plt.Line2D((), (), color="#38bdf8", marker="^", label="Lscale_up"),
            plt.Line2D((), (), color="#fb7185", marker="v", label="Lscale_down"),
        ],
        frameon=False,
        loc="upper left",
    )
    figure.savefig(path, dpi=180, facecolor="#0d1525")
    plt.close(figure)


def _write_q_figure(local, lower, upper, coherent, path: Path):
    q_vectors = {
        "lower source": standardized_displacement(lower),
        "local": standardized_displacement(local),
        "upper source": standardized_displacement(upper),
        "directional prior": standardized_displacement(coherent),
    }
    figure = plt.figure(figsize=(9.4, 7.0), constrained_layout=True)
    axis = figure.add_subplot(projection="3d")
    colors = ("#38bdf8", "#fbbf24", "#fb7185", "#a78bfa")
    for (label, vector), color in zip(q_vectors.items(), colors):
        axis.quiver(
            0, 0, 0, *vector, color=color, linewidth=2.4,
            arrow_length_ratio=0.08, label=label,
        )
    # q is dimensionless and the diagnosed values here are deliberately small.
    # Fit the axes to this comparison rather than imposing a unit-cube view
    # that makes the directional adjustment look like a point at the origin.
    extent = max(0.02, max(float(np.linalg.norm(vector)) for vector in q_vectors.values()) * 1.35)
    axis.set(xlim=(-extent, extent), ylim=(-extent, extent), zlim=(-extent, extent))
    axis.set_xlabel("q_w")
    axis.set_ylabel("q_rₜ")
    axis.set_zlabel("q_θₗ")
    axis.set_title("Dimensionless PDF-10 center-separation geometry")
    axis.view_init(elev=24, azim=-50)
    axis.legend(loc="upper left", frameon=False)
    figure.savefig(path, dpi=180, facecolor="#0d1525")
    plt.close(figure)


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    definition = base.CASE_DEFINITIONS["arm"]
    configuration = base._slider_configuration("arm")
    time_index = int(np.argmin(np.abs(configuration["times"] - TIME_SECONDS)))
    height_index = int(np.argmin(np.abs(configuration["heights"] - TARGET_HEIGHT_M)))
    time_seconds = int(configuration["times"][time_index])
    heights = configuration["heights"]
    tuning = base._tuning_from_values(tuple(control[5] for control in base.TUNING_CONTROLS))
    snapshot = base._raw_snapshot("arm", time_seconds, float(heights[height_index]))
    lower_snapshot = base._raw_snapshot("arm", time_seconds, float(heights[height_index - 1]))
    upper_snapshot = base._raw_snapshot("arm", time_seconds, float(heights[height_index + 1]))

    def diagnose(item):
        third = np.asarray((item.wp3, item.rtp3, item.thlp3), float)
        return diagnose_transport_2g_from_moments(item.mean, item.covariance, third, tuning)

    local, lower, upper = map(diagnose, (snapshot, lower_snapshot, upper_snapshot))
    lscale = compute_sam_lscale_profile(str(definition["run_dir"]), time_seconds)
    coherence = apply_local_vertical_coherence(
        snapshot.mean,
        snapshot.covariance,
        np.asarray((snapshot.wp3, snapshot.rtp3, snapshot.thlp3), float),
        tuning,
        local,
        lower=lower,
        lower_distance_m=snapshot.height_m - lower_snapshot.height_m,
        lower_reach_m=float(lscale.lscale_up_m[height_index - 1]),
        upper=upper,
        upper_distance_m=upper_snapshot.height_m - snapshot.height_m,
        upper_reach_m=float(lscale.lscale_down_m[height_index + 1]),
        settings=VerticalCoherenceSettings(enabled=True, max_blend=0.15),
    )

    rt_thl_path = FIGURE_DIR / "rt-thl-local-and-directional.png"
    lscale_path = FIGURE_DIR / "directional-lscale-reach.png"
    q_path = FIGURE_DIR / "q-coordinate-prior.png"
    _write_rt_thl_figure(snapshot, local, coherence.result, rt_thl_path)
    _write_lscale_figure(lscale, snapshot.height_m, lscale_path)
    _write_q_figure(local, lower, upper, coherence.result, q_path)

    revision = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    report = ReportBuilder(
        REPORT_ID,
        "Directional Lscale priors for PDF-10 vertical coherence",
        summary="Use CLUBB's parcel-derived Lscale_up/down as directional reach weights for a bounded prior on PDF-10's dimensionless center geometry, while retaining the selected level's local moments as the closure authority.",
        tags=("PDF-10", "Lscale", "SAM", "vertical coherence", "q coordinates"),
        replace=True,
        source_revision=revision,
    )
    report.heading("Recommendation")
    report.callout(
        "Adopt a directional one-interface prior first",
        "Restore the parcel-based mixing-length calculation as the default. For PDF-10, let a lower diagnosed PDF influence the next level only through the lower level's Lscale_up, and let an upper PDF influence it only through the upper level's Lscale_down. Blend that prior weakly into the local q geometry, then rerun the existing local PDF-10 reconstruction and PSD cap. Do not add a vertical matrix solve yet.",
        tone="success",
    )
    report.metrics(
        [
            ("Raw plane", f"ARM · {snapshot.elapsed_minutes:.0f} min · {snapshot.height_m:.0f} m", "selected reproducible raw-SAM snapshot"),
            ("Direct Lscale", f"{lscale.lscale_m[height_index]:.0f} m", "native calc_Lscale_directly at selected level"),
            ("Source reach", f"{lscale.lscale_up_m[height_index-1]:.0f} / {lscale.lscale_down_m[height_index+1]:.0f} m", "lower up / upper down"),
            ("Applied blend", f"{coherence.blend:.1%}", f"lower/upper support {coherence.lower_weight:.2f}/{coherence.upper_weight:.2f}"),
        ]
    )
    report.heading("Why q coordinates are useful")
    report.paragraph(
        "PDF-10 should not treat the two component centers as six unrelated physical numbers. Its local geometry is represented by one dimensionless separation vector q. With D equal to the diagonal matrix of the supplied grid standard deviations and mixture weights a and b, q generates both means while preserving the supplied grid mean exactly. This gives a coordinate system in which a 0.3-standard-deviation moisture separation and a 0.3-standard-deviation vertical-velocity separation can be compared, blended, and bounded on the same footing."
    )
    report.equation(r"\boldsymbol{\mu}_1=\overline{\mathbf{x}}+b\mathbf{D}\mathbf{q},\qquad \boldsymbol{\mu}_2=\overline{\mathbf{x}}-a\mathbf{D}\mathbf{q}", caption="The local grid mean is conserved by construction.")
    report.paragraph(
        "The between-component covariance is ab D q q^T D. The remaining covariance stays inside the component covariances, so changing q does not grant arbitrary physical freedom: the existing covariance budget and positive-semidefinite reconstruction still decide what is realizable. q is therefore the right object for a vertical prior: it describes relative PDF geometry, not an unmodified parcel anomaly being transported through changing background variance."
    )
    report.equation(r"\mathbf{C}=a\mathbf{\Sigma}_1+b\mathbf{\Sigma}_2+ab\,\mathbf{D}\mathbf{q}\mathbf{q}^{T}\mathbf{D}", caption="The same local covariance budget remains authoritative after the q adjustment.")
    report.figure(q_path, caption="At the selected ARM plane, immediate lower/upper q vectors supply a weak directional geometric prior. The reconstructed PDF remains tied to the target level's own mean and covariance.")
    report.heading("Directional Lscale prior")
    report.paragraph(
        "The native direct mixing-length calculation follows buoyant, entraining parcel trajectories. Lscale_up is the upward distance over which a parcel launched at a lower level exhausts its turbulent-energy budget; Lscale_down is the analogous downward distance from an upper source. Those are naturally directional influence lengths, unlike the tau-derived scalar length that does not carry a source direction."
    )
    report.equation(r"W_{k-1\rightarrow k}=G\exp[-(\Delta z/L_{\mathrm{up},k-1})^2],\qquad W_{k+1\rightarrow k}=G\exp[-(\Delta z/L_{\mathrm{down},k+1})^2]", caption="A source-sensitive, one-interface reach weight. G is the existing q-direction/activity gate.")
    report.paragraph(
        "Combine the two weighted neighbor geometries to form q_prior, then use q_new = (1-alpha) q_local + alpha q_prior, with alpha capped at 0.25. The Notes implementation uses a 15% default cap. Finally, it calls the ordinary local PDF-10 diagnosis with q_new; widths, third-moment handling, covariance conservation, and PSD safety remain local algebra."
    )
    report.figure(lscale_path, caption="Actual native direct-Lscale arrows from the same raw-SAM column. Horizontal coordinate is source altitude; vertical coordinate is the source plus/minus the directional reach. The gold line marks the selected plane.")
    report.heading("What changes in the selected r_t–theta_l plane")
    report.paragraph(
        "The raw r_t–theta_l heatmap remains the evidence. Dashed component loops are the independent PDF-10 diagnosis; solid loops are the same local reconstruction after the directional q prior. A useful change should make component geometry more vertically consistent without moving the distribution outside the raw thermodynamic shape or consuming unsupported covariance."
    )
    report.figure(rt_thl_path, caption="Raw ARM r_t–theta_l probability with local versus directional-prior component footprints. The directional step changes only q before the same selected-plane PDF-10 reconstruction.")
    report.heading("Inputs, limits, and next test")
    report.paragraph(
        "This SAM wrapper calls the compiled Python API's calc_Lscale_directly; it does not reproduce parcel integration in Python. It precomputes one full-column resolved profile at each selected SAM time: means of r_t, theta_l, and r_c; r_t/theta_l covariance; pressure; and w variance. The micro files contain no u or v, so the one explicit approximation is em = 1.5 w'^2. This is suitable for comparing directional structure, not for claiming an exact host-model TKE diagnosis."
    )
    report.callout(
        "Next evaluation",
        "Run local-off versus directional-on placement scores over ARM low-level V regimes and elevated negative-correlation regimes. Score raw-SAM center/footprint continuity and cloud/transport diagnostics separately. Keep the direct Lscale path if it improves the former without degrading the latter. Only if one-interface priors are demonstrably too local should we test a banded vertical precision solve; that would be a second phase, not a prerequisite.",
        tone="note",
    )
    report.code_file(Path(__file__), language="python", caption="Reproduction script; it uses the native F2PY calc_Lscale_directly wrapper and embeds the generated PNGs.")
    data = {
        "source_file": snapshot.source_file,
        "time_seconds": snapshot.elapsed_seconds,
        "height_m": snapshot.height_m,
        "lscale_m": float(lscale.lscale_m[height_index]),
        "lower_source_lscale_up_m": float(lscale.lscale_up_m[height_index - 1]),
        "upper_source_lscale_down_m": float(lscale.lscale_down_m[height_index + 1]),
        "blend": coherence.blend,
        "lower_weight": coherence.lower_weight,
        "upper_weight": coherence.upper_weight,
        "approximation": lscale.approximation,
    }
    (report.stage / "data").mkdir(exist_ok=True)
    (report.stage / "data" / "directional-lscale-snapshot.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8"
    )
    report.publish()


if __name__ == "__main__":
    main()
