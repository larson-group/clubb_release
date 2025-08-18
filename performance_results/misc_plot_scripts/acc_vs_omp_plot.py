import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import kaleido

# Define input files per device and directive method
files = {
    "V100 OpenACC": "timing_results_gptl_new/V100_nvhpc_noasync_1x1_casper_arm.csv",
    "V100 OpenMP":  "timing_results_gptl_new/V100_nvhpc_omp_1x1_casper_arm.csv",
    "A100 OpenACC": "timing_results_gptl_new/A100_nvhpc_noasnyc_1x1_derecho_arm.csv",
    "A100 OpenMP":  "timing_results_gptl_new/A100_nvhpc_omp_1x1_derecho_arm.csv",
    "H100 OpenACC": "timing_results_gptl_new/H100_nvhpc_noasync_1x1_casper_arm.csv",
    "H100 OpenMP":  "timing_results_gptl_new/H100_nvhpc_omp_1x1_casper_arm.csv"
}

# Marker symbols for visual distinction
symbol_sequence = ["circle", "square", "diamond"]

# Load data
data = {}
for label, path in files.items():
    df = pd.read_csv(path, comment="#")
    df = df[df["ngrdcol"] >= 4]
    df.sort_values("ngrdcol", inplace=True)
    data[label] = df[["ngrdcol", "compute_i"]].copy()

# Plot setup
fig = go.Figure()
devices = ["V100", "A100", "H100"]

for idx, device in enumerate(devices):
    acc_label = f"{device} OpenACC"
    omp_label = f"{device} OpenMP"

    if acc_label not in data or omp_label not in data:
        continue

    merged = pd.merge(
        data[acc_label],
        data[omp_label],
        on="ngrdcol",
        suffixes=("_acc", "_omp")
    )
    merged["speedup"] = merged["compute_i_acc"] / merged["compute_i_omp"]

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
    annotation_text="Baseline (OpenACC)",
    annotation_position="top right"
)

font_size = 16

# Layout configuration
fig.update_layout(
    title="Directive Method: OpenMP vs OpenACC",
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
        title="Throughput Ratio (OpenMP / OpenACC)",
        title_font=dict(size=font_size),
        range=[0.4, 1.4],
        rangemode="normal",
        linecolor='black',
        linewidth=1,
        mirror=True
    ),
    legend=dict(
        x=0.01, y=0.01,
        xanchor="left", yanchor="bottom",
        font=dict(size=font_size),
        bgcolor="rgba(255,255,255,0.5)"
    ),
    width=1600,
    height=800,
    margin=dict(l=10, r=10, t=50, b=10),
    template="plotly_white"
)

# Export as PNG
pio.write_image(
    fig,
    "openacc_vs_openmp_speedup.png",
    width=800,
    height=400,
    scale=3.125
)

#fig.show()
