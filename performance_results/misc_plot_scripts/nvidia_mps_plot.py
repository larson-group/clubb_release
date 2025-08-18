import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import kaleido

# File paths for different configurations
files = {
    "A100_4x4":  "timing_results_gptl_new/A100_nvhpc_noasync_4x4_derecho_arm.csv",
    "A100_4x8":  "timing_results_gptl_new/A100_nvhpc_noasync_4x8_134nz_derecho_arm.csv",
    "A100_4x16": "timing_results_gptl_new/A100_nvhpc_noasync_4x16_134nz_derecho_arm.csv",
    "A100_4x32": "timing_results_gptl_new/A100_nvhpc_noasync_4x32_134nz_derecho_arm.csv",
    "A100_4x64": "timing_results_gptl_new/A100_nvhpc_noasync_4x64_134nz_derecho_arm.csv"
}

# Marker symbol sequence
symbol_sequence = ["circle", "square", "diamond", "cross", "x", "triangle-up", "star"]

# Load and prepare data
data = {}
for label, path in files.items():
    df = pd.read_csv(path, comment="#")
    df.sort_values("ngrdcol", inplace=True)
    data[label] = df[["ngrdcol", "compute_i"]].copy()

# Use 4x4 as baseline
baseline_label = "A100_4x4"
baseline_df = data[baseline_label].rename(columns={"compute_i": "compute_i_baseline"})

# Create figure
fig = go.Figure()

# Add traces with distinct marker symbols
symbol_index = 0
for label, df in data.items():
    if label == baseline_label:
        continue  # skip baseline trace

    # Merge with baseline on ngrdcol
    merged = pd.merge(baseline_df, df, on="ngrdcol")
    merged["speedup"] = merged["compute_i_baseline"] / merged["compute_i"]

    fig.add_trace(go.Scatter(
        x=merged["ngrdcol"],
        y=merged["speedup"],
        mode="lines+markers",
        name=label,
        marker=dict(size=12, symbol=symbol_sequence[symbol_index % len(symbol_sequence)])
    ))

    symbol_index += 1

# Add baseline line at y=1
fig.add_hline(
    y=1,
    line_dash="dash",
    line_color="gray",
    annotation_text="Baseline (4x4)",
    annotation_position="bottom left"
)

font_size = 16

# Layout
fig.update_layout(
    title="Oversubscribing GPUs with Nvidia MPS",
    title_font=dict(size=font_size + 4),
    font=dict(size=font_size),
    xaxis=dict(
        title="Batch Size (columns)",
        title_font=dict(size=font_size),
        type="log",
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    yaxis=dict(
        title="Throughput Ratio (Oversubscribed / 1 Task per GPU)",
        title_font=dict(size=font_size-2),
        range=[0.8, 1.35],
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    legend=dict(
        x=0.99, y=0.01,
        xanchor="right", yanchor="bottom",
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
    "nvidia_mps_speedup.png",
    width=800,
    height=400,
    scale=3.125
)

#fig.show()
