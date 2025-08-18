import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import kaleido

# File paths for async and non-async runs per GPU
files = {
    "Async NVHPC V100": "timing_results_gptl_new/V100_nvhpc_async_1x1_casper_arm.csv",
    "Non‑Async V100":   "timing_results_gptl_new/V100_nvhpc_noasync_1x1_casper_arm.csv",
    "Async NVHPC A100": "timing_results_gptl_new/A100_nvhpc_asnyc_1x1_derecho_arm.csv",
    "Non‑Async A100":   "timing_results_gptl_new/A100_nvhpc_noasnyc_1x1_derecho_arm.csv",
    "Async NVHPC H100": "timing_results_gptl_new/H100_nvhpc_async_1x1_derecho_arm.csv",
    "Non‑Async H100":   "timing_results_gptl_new/H100_nvhpc_noasync_1x1_casper_arm.csv"
}

# Marker symbols for each GPU
symbol_sequence = ["circle", "square", "diamond"]

# Load and prepare data
data = {}
for label, path in files.items():
    df = pd.read_csv(path, comment="#")
    df.sort_values("ngrdcol", inplace=True)
    df = df[df["ngrdcol"] >= 4]  # keep only rows where ngrdcol >= 4
    data[label] = df[["ngrdcol", "compute_i"]].copy()

# Merge async and non-async results per device and calculate speedup
fig = go.Figure()
devices = ["V100", "A100", "H100"]

for idx, device in enumerate(devices):
    async_label = f"Async NVHPC {device}"
    nonasync_label = f"Non‑Async {device}"

    merged = pd.merge(
        data[async_label],
        data[nonasync_label],
        on="ngrdcol",
        suffixes=("_async", "_base")
    )
    merged["speedup"] = merged["compute_i_base"] / merged["compute_i_async"]

    fig.add_trace(go.Scatter(
        x=merged["ngrdcol"],
        y=merged["speedup"],
        mode="lines+markers",
        name=device+"_1x1",
        marker=dict(symbol=symbol_sequence[idx % len(symbol_sequence)], size=14)
    ))

# Add baseline line at y=1
fig.add_hline(
    y=1,
    line_dash="dash",
    line_color="gray",
    annotation_text="Baseline (Non-async)",
    annotation_position="top left"
)

font_size = 16

# Update layout
fig.update_layout(
    title="Async vs Non-Async Throughput Ratio",
    title_font=dict(size=font_size+4),
    font=dict(size=font_size),
    xaxis=dict(
        title="Batch Size",
        title_font=dict(size=font_size),
        type="log",
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    yaxis=dict(
        title="Throughput Ratio (Async / Non‑Async)",
        title_font=dict(size=font_size),
        range=[0.8, 2.1],
        rangemode="normal",
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    legend=dict(
        x=0.99, y=0.99,
        xanchor="right", yanchor="top",
        font=dict(size=font_size),
        bgcolor="rgba(255,255,255,0.5)"
    ),
    width=1600,
    height=800,
    margin=dict(l=10, r=10, t=50, b=10),
    template="plotly_white"
)

# Export to PNG
pio.write_image(
    fig,
    "async_speedup_comparison.png",
    width=800,
    height=400,
    scale=3.125
)

#fig.show()
