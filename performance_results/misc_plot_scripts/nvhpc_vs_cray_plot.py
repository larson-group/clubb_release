import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import kaleido

# Load data
files = {
    "A100 NVHPC": "timing_results_gptl_new/A100_nvhpc_noasync_4x4_derecho_arm.csv",
    "A100 CRAY":  "timing_results_gptl_new/A100_cray_noasync_4x4_derecho_arm.csv"
}

data = {}
for label, path in files.items():
    df = pd.read_csv(path, comment="#")
    df.sort_values("ngrdcol", inplace=True)
    data[label] = df[["ngrdcol", "compute_i"]].copy()

# Merge and compute speedup (NVHPC as baseline)
merged = pd.merge(
    data["A100 NVHPC"],
    data["A100 CRAY"],
    on="ngrdcol",
    suffixes=("_nvhpc", "_cray")
)
merged["speedup"] = merged["compute_i_nvhpc"] / merged["compute_i_cray"]

# Determine y-axis buffer
ymin = merged["speedup"].min()
y_buffer = max((1 - ymin) * 0.1, 0.05)
y_lower = min(ymin, 1 - y_buffer)

# Create figure
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=merged["ngrdcol"],
    y=merged["speedup"],
    mode="lines+markers",
    name="Speedup",
    marker=dict(symbol="circle", size=14)
))

# Add a horizontal line at y=1
fig.add_hline(
    y=1,
    line_dash="dash",
    line_color="gray",
    annotation_text="Baseline (NVHPC)",
    annotation_position="top left"
)

font_size = 16

# Update layout
fig.update_layout(
    title="Compiler Comparison: NVHPC vs Cray",
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
        title="Throughput Ratio (CRAY / NVHPC)",
        title_font=dict(size=font_size),
        range=[0.8, 1.4],
        rangemode="normal",
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    legend=dict(
        x=0.01, y=0.99,
        xanchor="left", yanchor="top",
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
    "nvhpc_vs_cray_speedup.png",
    width=800,
    height=400,
    scale=3.125
)

#fig.show()
