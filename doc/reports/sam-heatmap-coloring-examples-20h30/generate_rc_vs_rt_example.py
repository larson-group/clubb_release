#!/usr/bin/env python3
"""Render five matched raw-SAM w-r_t color-signal panels for this report.

The saved PNG is a report artifact, not a runtime dashboard product.  Re-run
this script only when deliberately updating the example, then commit both the
script and its PNG/JSON outputs together with the report HTML that embeds them.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import numpy as np


REPORT_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPORT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_dash_report_matplotlib")

from dash_app.shared.bivariate_heatmap import (
    BIVARIATE_LEVELS,
    robust_field_upper,
    signed_bivariate_codes,
    signed_bivariate_reference_rgb,
)
from utilities.sam_3d_reference import DEFAULT_SAM_RUN, load_snapshot


ELAPSED_SECONDS = 294 * 60
TARGET_HEIGHT_M = 820.0
BINS = 72


def _axis_edges(rt_g_per_kg: np.ndarray, w: np.ndarray, rc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Use one physical w-r_t grid for both color signals."""
    x_low, x_high = np.percentile(rt_g_per_kg, (0.1, 99.9))
    y_low, y_high = np.percentile(w, (0.1, 99.9))
    cloudy = rc > 0.0
    if np.any(cloudy):
        x_low = min(x_low, float(np.percentile(rt_g_per_kg[cloudy], 0.5)))
        x_high = max(x_high, float(np.percentile(rt_g_per_kg[cloudy], 99.5)))
        y_low = min(y_low, float(np.percentile(w[cloudy], 0.5)))
        y_high = max(y_high, float(np.percentile(w[cloudy], 99.5)))
    x_pad = max(0.04 * (x_high - x_low), 1.0e-3)
    y_pad = max(0.04 * (y_high - y_low), 1.0e-3)
    return (
        np.linspace(x_low - x_pad, x_high + x_pad, BINS + 1),
        np.linspace(y_low - y_pad, y_high + y_pad, BINS + 1),
    )


def _binned_field(x, w, x_edges, y_edges, local_signal):
    probability = np.histogram2d(x, w, bins=(x_edges, y_edges))[0].T / x.size
    signal = np.histogram2d(x, w, bins=(x_edges, y_edges), weights=local_signal)[0].T / x.size
    return probability, signal


def _rgb(probability, signed_signal, probability_upper, signal_upper):
    codes, _, _ = signed_bivariate_codes(
        probability, signed_signal, probability_upper, signal_upper
    )
    signed_levels = 2 * BIVARIATE_LEVELS - 1
    palette = np.empty((BIVARIATE_LEVELS * signed_levels, 3), dtype=np.uint8)
    for probability_bin in range(BIVARIATE_LEVELS):
        for signed_bin in range(signed_levels):
            code = probability_bin * signed_levels + signed_bin
            palette[code] = signed_bivariate_reference_rgb(
                probability_bin / (BIVARIATE_LEVELS - 1),
                (signed_bin - (BIVARIATE_LEVELS - 1)) / (BIVARIATE_LEVELS - 1),
            )
    return palette[codes]


def main() -> int:
    from matplotlib import pyplot as plt

    snapshot = load_snapshot(DEFAULT_SAM_RUN, ELAPSED_SECONDS, TARGET_HEIGHT_M)
    w = np.asarray(snapshot.samples[:, 0], dtype=float)
    rt = np.asarray(snapshot.samples[:, 1], dtype=float) * 1000.0
    rc = np.asarray(snapshot.rc_samples, dtype=float) * 1000.0
    chi = np.asarray(snapshot.chi_samples, dtype=float) * 1000.0
    thl = np.asarray(snapshot.samples[:, 2], dtype=float)
    x_edges, y_edges = _axis_edges(rt, w, rc)

    w_prime = w - np.mean(w)
    rc_prime = rc - np.mean(rc)
    rt_prime = rt - np.mean(rt)
    chi_prime = chi - np.mean(chi)
    thl_prime = thl - np.mean(thl)
    probability, signal_rc = _binned_field(
        rt, w, x_edges, y_edges, w_prime * rc_prime
    )
    signal_specs = (
        (
            "wprcp",
            signal_rc,
            "Current cloud-water covariance",
            "w′rᶜ′  [m s⁻¹ g kg⁻¹]",
        ),
        (
            "wprtp",
            _binned_field(rt, w, x_edges, y_edges, w_prime * rt_prime)[1],
            "Total-water covariance",
            "w′rₜ′  [m s⁻¹ g kg⁻¹]",
        ),
        (
            "wpchi",
            _binned_field(rt, w, x_edges, y_edges, w_prime * chi_prime)[1],
            "Saturation-excess covariance",
            "w′χ′  [m s⁻¹ g kg⁻¹]",
        ),
        (
            "wpthlp",
            _binned_field(rt, w, x_edges, y_edges, w_prime * thl_prime)[1],
            "Liquid-water-potential-temperature covariance",
            "w′θₗ′  [m K s⁻¹]",
        ),
        (
            "wprtp2",
            _binned_field(rt, w, x_edges, y_edges, w_prime * np.square(rt_prime))[1],
            "Total-water third-moment contribution",
            "w′rₜ′²  [m s⁻¹ (g kg⁻¹)²]",
        ),
    )
    probability_upper = robust_field_upper((probability,))
    signal_uppers = {
        name: robust_field_upper((np.abs(signal),))
        for name, signal, _title, _formula in signal_specs
    }

    figure, axes = plt.subplots(2, 3, figsize=(18.0, 10.4), dpi=160, facecolor="#0d1525")
    for index, (name, signal, title, formula) in enumerate(signal_specs):
        axis = axes.flat[index]
        upper = signal_uppers[name]
        axis.set_facecolor("#0d1525")
        axis.imshow(
            _rgb(probability, signal, probability_upper, upper),
            origin="lower",
            extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
            aspect="auto",
            interpolation="nearest",
        )
        axis.axhline(0.0, color="#cbd5e1", alpha=0.38, linewidth=0.9, linestyle="--")
        axis.set_title(f"{title}\n{formula}", color="#e8edf5", fontsize=11, pad=11)
        axis.set_xlabel("Total water, rₜ [g/kg]", color="#e8edf5")
        if index % 3 == 0:
            axis.set_ylabel("Vertical velocity, w [m/s]", color="#e8edf5")
        axis.tick_params(colors="#cbd5e1")
        for spine in axis.spines.values():
            spine.set_color("#40516a")
    axes.flat[-1].set_facecolor("#0d1525")
    axes.flat[-1].axis("off")
    axes.flat[-1].text(
        0.08,
        0.16,
        "One fixed ARM raw-SAM plane",
        color="#f8fafc",
        fontsize=13,
        fontweight="bold",
        transform=axes.flat[-1].transAxes,
    )
    axes.flat[-1].text(
        0.08,
        0.92,
        "Gold = probability\nBlue/red = positive/negative\nlocal moment contribution\n\nEach panel has its own robust\n99.5th-percentile magnitude cap.",
        color="#cbd5e1",
        fontsize=10,
        linespacing=1.55,
        transform=axes.flat[-1].transAxes,
        va="top",
    )
    figure.suptitle(
        f"Raw ARM SAM: {snapshot.elapsed_minutes:.0f} min / {snapshot.height_m:.0f} m — same w–rₜ bins and probability, five color signals",
        color="#f8fafc",
        fontsize=14,
        y=0.98,
    )
    figure.text(
        0.5,
        0.012,
        "All covariance panels center w and q.  The final panel is the third-order kernel w′rₜ′², so it is not a covariance.",
        ha="center",
        color="#a9b8cd",
        fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.04, 1, 0.93))
    figure_dir = REPORT_DIR / "figures"
    figure_dir.mkdir(exist_ok=True)
    figure.savefig(figure_dir / "raw-sam-wrt-color-signals.png", facecolor=figure.get_facecolor())
    plt.close(figure)

    metadata = {
        "source_file": snapshot.source_file,
        "elapsed_minutes": snapshot.elapsed_minutes,
        "height_m": snapshot.height_m,
        "sample_count": snapshot.sample_count,
        "bins": BINS,
        "probability_upper": probability_upper,
        "signals": {
            name: {
                "formula": formula,
                "robust_positive_upper": signal_uppers[name],
                "integrated_moment": float(np.sum(signal)),
            }
            for name, signal, _title, formula in signal_specs
        },
    }
    (REPORT_DIR / "data").mkdir(exist_ok=True)
    (REPORT_DIR / "data" / "raw-sam-wrt-color-signals.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
