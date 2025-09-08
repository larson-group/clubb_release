
import os
import re
import glob
import pandas as pd
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
import numpy as np

# Step 1: Aggregate all metric files
metric_files = sorted(
    glob.glob("Intel6430_uarch_vtune_multicore/vtune_*col/metrics_*col.csv"),
    key=lambda x: int(re.search(r'vtune_(\d+)col', x).group(1))
)

# Step 2: Parse and sum metrics per file
totals_by_ncol = {}
for file_path in metric_files:
    match = re.search(r'vtune_(\d+)col', file_path)
    if not match:
        continue
    ncol = int(match.group(1))
    df = pd.read_csv(file_path, sep='\t', engine='python')
    metric_cols = [col for col in df.columns if col.startswith("Hardware Event Count:")]
    cleaned_cols = [col.replace("Hardware Event Count:", "").strip() for col in metric_cols]
    summed_metrics = df[metric_cols].sum()
    summed_metrics.index = cleaned_cols
    totals_by_ncol[ncol] = summed_metrics

# Step 3: Create combined DataFrame
totals_df = pd.DataFrame.from_dict(totals_by_ncol, orient='index').sort_index()
totals_df.index.name = 'Ncol'

# Step 3.1: Compute derived metrics
totals_df['CPI'] = totals_df['CPU_CLK_UNHALTED.THREAD'] / totals_df['INST_RETIRED.ANY']
totals_df['IPC'] = totals_df['INST_RETIRED.ANY'] / totals_df['CPU_CLK_UNHALTED.THREAD']
totals_df['Vectorization Ratio'] = totals_df['FP_ARITH_INST_RETIRED.VECTOR'] / (
    totals_df['FP_ARITH_INST_RETIRED.VECTOR'] + totals_df['FP_ARITH_INST_RETIRED.SCALAR']
)

# --- Add after Step 3.1 (after CPI/Vectorization Ratio) ---
stall_raw_cols = [
    "MEMORY_ACTIVITY.STALLS_L1D_MISS",
    "MEMORY_ACTIVITY.STALLS_L2_MISS",
    "MEMORY_ACTIVITY.STALLS_L3_MISS",
]

# Ensure columns exist even if some are absent in certain files
for col in stall_raw_cols:
    if col not in totals_df.columns:
        totals_df[col] = 0

# Per-cycle normalization (percentage of total core cycles)
denom = totals_df["CPU_CLK_UNHALTED.THREAD"].replace(0, np.nan)  # avoid divide-by-zero
totals_df["Stall L1D Miss %"] = 100.0 * (totals_df["MEMORY_ACTIVITY.STALLS_L1D_MISS"] / denom)
totals_df["Stall L2 Miss %"] = 100.0 * (totals_df["MEMORY_ACTIVITY.STALLS_L2_MISS"] / denom)
totals_df["Stall L3 Miss %"] = 100.0 * (totals_df["MEMORY_ACTIVITY.STALLS_L3_MISS"] / denom)

# If denom was NaN (no cycles), fill NaNs introduced by division with 0 for plotting
totals_df[["Stall L1D Miss %","Stall L2 Miss %","Stall L3 Miss %"]] = \
    totals_df[["Stall L1D Miss %","Stall L2 Miss %","Stall L3 Miss %"]].fillna(0.0)



graph_config = {
    'toImageButtonOptions': {
        'format': 'png',
        'filename': 'custom_image',
        'height': 400,
        'width': 800,
        'scale': 3.125  # 300 DPI equivalent
    }
}


# Step 4: Build Dash App
app = dash.Dash(__name__)
app.layout = html.Div([
    html.H1("VTune Metrics Dashboard"),

    html.Div([
        # Left column: two stacked graphs
        html.Div([
            dcc.Graph(id='metric-graph', config=graph_config),
            dcc.Graph(
                id='vectorization-graph',
                figure={
                    'data': [
                        go.Scatter(
                            x=totals_df.index,
                            y=totals_df['Vectorization Ratio'],
                            mode='lines+markers',
                            name='Vectorization Ratio',
                            marker=dict(size=6)  
                        )
                    ],
                    'layout': go.Layout(
                        title=dict(
                            text="CPU vectorization ratio",
                            font=dict(size=24),  # Title font size
                            x=0.5,  # Center horizontally
                            xanchor="center"  # Anchor relative to the center
                        ),
                        xaxis=dict(
                            title=dict(
                                text="Per Core Batch Size (columns per core)",
                                font=dict(size=18)
                            ),
                            showline=True,
                            linecolor="black",
                            linewidth=1,
                            mirror=True,
                            dtick=4
                        ),
                        yaxis=dict(
                            title=dict(
                                text="Instruction Ratio (vector FLOPs / total FLOPs)",
                                font=dict(size=16)
                            ),
                            showline=True,
                            linecolor="black",
                            linewidth=1,
                            mirror=True,
                            dtick=0.05
                        ),
                        margin=dict(l=10, r=10, t=50, b=10),
                        template="plotly_white"
                    )
                },
                config=graph_config
            ),
            dcc.Graph(
                id='stall-graph',
                figure={
                    'data': [
                        go.Scatter(
                            x=totals_df.index,
                            y=totals_df['Stall L3 Miss %'],
                            name='STALLS_L3_MISS',
                            mode='lines',
                            line=dict(width=0.5, color="red"),
                            stackgroup='one',       # enables stacking
                            groupnorm=None,         # or 'percent' if you want percentages
                        ),
                        go.Scatter(
                            x=totals_df.index,
                            y=totals_df['Stall L2 Miss %'],
                            name='STALLS_L2_MISS',
                            mode='lines',
                            line=dict(width=0.5, color="yellow"),
                            stackgroup='one',
                        ),
                        go.Scatter(
                            x=totals_df.index,
                            y=totals_df['Stall L1D Miss %'],
                            name='STALLS_L1D_MISS',
                            mode='lines',
                            line=dict(width=0.5, color="green"),
                            stackgroup='one',
                        ),
                    ],
                    'layout': go.Layout(
                        title=dict(
                            text="CPU stalls due to cache misses",
                            #text="Cache miss stalls on Intel 6430 CPU",
                            font=dict(size=24),  # Title font size
                            x=0.5,  # Center horizontally
                            xanchor="center"  # Anchor relative to the center
                        ),
                        xaxis=dict(
                            title=dict(
                                text="Per Core Batch Size (columns per core)",
                                font=dict(size=18)
                            ),
                            showline=True,
                            linecolor="black",
                            linewidth=1,
                            mirror=True,
                            dtick=4
                        ),
                        yaxis=dict(
                            title=dict(
                                text="Stalled Cycles (% of total cycles)",
                                font=dict(size=18)
                            ),
                            showline=True,
                            linecolor="black",
                            linewidth=1,
                            mirror=True
                        ),
                        margin=dict(l=10, r=10, t=50, b=10),
                        template="plotly_white",
                        legend=dict(
                            orientation="v",
                            yanchor="top",
                            y=0.98,
                            xanchor="left",
                            x=0.02,
                            bgcolor="rgba(255,255,255,0.5)"
                        )
                    )
                },
                config=graph_config
            ),

        ], style={
            'width': '74%',
            'display': 'inline-block',
            'verticalAlign': 'top'
        }),

        # Right column: checklist
        html.Div([
            html.Label("Select Metrics:"),
            dcc.Checklist(
                id='metric-selector',
                options=[{'label': col, 'value': col} 
                        for col in ['CPI', 'IPC'] + 
                        [col for col in totals_df.columns if col not in ['CPI', 'IPC', 'Vectorization Ratio']]],
                value=['CPI'],
                style={'height': '90vh', 'overflowY': 'scroll'}
            )
        ], style={
            'width': '24%',
            'paddingLeft': '2%',
            'display': 'inline-block',
            'verticalAlign': 'top'
        })

    ], style={
        'display': 'flex',
        'flexDirection': 'row',
        'alignItems': 'flex-start'
    })
])

@app.callback(
    Output('metric-graph', 'figure'),
    [Input('metric-selector', 'value')]
)
def update_graph(selected_metrics):
    fig = go.Figure()

    # Plot selected metrics
    for metric in selected_metrics:
        fig.add_trace(go.Scatter(
            x=totals_df.index,
            y=totals_df[metric],
            mode='lines+markers',
            marker=dict(size=8),
            name=metric
        ))

    # Determine which special metrics are present
    has_cpi = "CPI" in selected_metrics
    has_ipc = "IPC" in selected_metrics

    # Add reference lines (once)
    if has_cpi:
        fig.add_hline(
            y=0.25,
            line=dict(color="grey", dash="dash"),
            annotation_text="CPI lower bound (0.25)",
            annotation_position="top right"
        )
    if has_ipc:
        fig.add_hline(
            y=4.0,
            line=dict(color="grey", dash="dash"),
            annotation_text="IPC upper bound (4)",
            annotation_position="bottom right"
        )

    # Explicit y-axis bounds
    y_range = None
    if has_cpi and has_ipc:
        y_range = [0, 6]       # covers both
    elif has_cpi:
        y_range = [0, 6]
    elif has_ipc:
        y_range = [0, 4.5]

    # Conditional y-axis title
    if has_cpi and has_ipc:
        y_title = "Cycles per Instruction (CPI) and Instructions per Cycle (IPC)"
    elif has_cpi:
        y_title = "Cycles per Instruction (CPI)"
    elif has_ipc:
        y_title = "Instructions per Cycle (IPC)"
    else:
        y_title = "Metric Value"

    fig.update_layout(
        title=dict(
            text="CLUBB VTune Metric on Intel 6430 CPU",
            font=dict(size=24)
        ),
        xaxis=dict(
            title=dict(text="Per Core Batch Size (columns per core)", font=dict(size=18)),
            showline=True, linecolor="black", linewidth=1, mirror=True
        ),
        yaxis=dict(
            title=dict(text=y_title, font=dict(size=18)),
            showline=True, linecolor="black", linewidth=1, mirror=True,
            range=y_range
        ),
        margin=dict(l=10, r=10, t=50, b=10),
        template="plotly_white"
    )
    return fig


if __name__ == '__main__':
    app.run(debug=True)
