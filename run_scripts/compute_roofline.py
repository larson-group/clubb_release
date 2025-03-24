import csv
import sys
import numpy as np
import subprocess
import pandas as pd
import dash
from dash import dcc, html
import plotly.graph_objs as go


# Verify command-line arguments
if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <ncu csv file>")
    sys.exit(1)

input_path = sys.argv[1]
output_path = "roofline_metrics.csv" #sys.argv[2]

# run_cmd = f"ncu --import {input_path} --csv --page raw > {output_path}"
# result = subprocess.run(run_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Prepare indices for needed columns
kernel_idx = time_idx = bytes_idx = sp_idx = dp_idx = hp_idx = None

# Read the input CSV (skip the units header) and identify column indices
with open(input_path, 'r') as infile:

    reader = csv.reader(infile)
    header_names = next(reader)        # first header line: metric names
    header_units = next(reader)        # second header line: metric units (can be ignored after index check)

    try:
        kernel_idx = header_names.index("Kernel Name")
        time_idx   = header_names.index("gpu__time_duration.sum")
        bytes_idx  = header_names.index("dram__bytes.sum.per_second")
        sp_idx     = header_names.index("derived__sm__sass_thread_inst_executed_op_ffma_pred_on_x2")
        dp_idx     = header_names.index("derived__sm__sass_thread_inst_executed_op_dfma_pred_on_x2")
        hp_idx     = header_names.index("derived__sm__sass_thread_inst_executed_op_hfma_pred_on_x2")
    except ValueError as e:
        sys.stderr.write(f"Error: Required column not found in CSV header ({e})\n")
        sys.exit(1)

    kernel_count = {}

    # Dictionary to accumulate metrics per kernel
    kernel_data = {}
    for row in reader:

        if not row or len(row) <= time_idx:
            continue  # skip empty or malformed lines

        #gpu_freq = int(float(row[header_names.index("device__attribute_max_gpu_frequency_khz")])) * 1e3
        cycles = int(float(row[header_names.index("gpc__cycles_elapsed.max")]))


        # Parse metrics
        time_ms = float(row[time_idx])
        bytes_per_sec = float(row[bytes_idx])  # in TBytes/s

        flops_sp = 0
        flops_sp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_fadd_pred_on.sum.per_cycle_elapsed")]))
        flops_sp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_fmul_pred_on.sum.per_cycle_elapsed")]))
        flops_sp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_ffma_pred_on.sum.per_cycle_elapsed")])) * 2
        flops_sp *= cycles


        flops_dp = 0
        flops_dp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_dadd_pred_on.sum.per_cycle_elapsed")]))
        flops_dp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_dmul_pred_on.sum.per_cycle_elapsed")]))
        flops_dp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_dfma_pred_on.sum.per_cycle_elapsed")])) * 2
        flops_dp *= cycles


        flops_hp = 0
        flops_hp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_hadd_pred_on.sum.per_cycle_elapsed")]))
        flops_hp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_hmul_pred_on.sum.per_cycle_elapsed")]))
        flops_hp += int(float(row[header_names.index("smsp__sass_thread_inst_executed_op_hfma_pred_on.sum.per_cycle_elapsed")])) * 2
        flops_hp *= cycles


        kname = row[kernel_idx]

        if kname in kernel_count:
          kernel_count[kname] += 1
        else:
          kernel_count[kname] = 1
        
        kname_unique = f"{kname}_{kernel_count[kname]}"

        # Convert time to seconds and compute total bytes for this kernel launch
        time_sec = time_ms * 1e-3
        total_bytes = bytes_per_sec * time_sec * 1e12  # TByte/s * seconds -> bytes

        # Initialize kernel entry if not present
        if kname_unique not in kernel_data:
            kernel_data[kname_unique] = {
                "count": 0,
                "total_time": 0.0,
                "total_bytes": 0.0,
                "total_sp": 0,
                "total_dp": 0,
                "total_hp": 0
            }

        entry = kernel_data[kname_unique]
        entry["count"] += 1
        entry["total_time"] += time_sec
        entry["total_bytes"] += total_bytes
        entry["total_sp"] += flops_sp
        entry["total_dp"] += flops_dp#+flops_sp+flops_hp
        entry["total_hp"] += flops_hp

        # entry = kernel_data[kname_unique]
        # entry["count"] = 1
        # entry["total_time"] += max(entry["total_time"], time_sec)
        # entry["total_bytes"] += max(entry["total_bytes"], total_bytes)
        # entry["total_sp"] += max(entry["total_sp"], flops_sp)
        # entry["total_dp"] += max(entry["total_dp"], flops_dp)#+flops_sp+flops_hp)
        # entry["total_hp"] += max(entry["total_hp"], flops_hp)

# Write the aggregated results to the output CSV
with open(output_path, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    # Output header
    writer.writerow([
        "Kernel Name",
        "Number of Calls",
        "Average Duration (s)",
        "Total Bytes",
        "Total FLOPs (SP)",
        "Total FLOPs (DP)",
        "Total FLOPs (HP)",
        "AI (SP)",
        "AI (DP)",
        "AI (HP)",
        "Achieved Performance (SP)",
        "Achieved Performance (DP)",
        "Achieved Performance (HP)"
    ])
    # Write one row per kernel with computed metrics
    for kname, agg in kernel_data.items():
        count = agg["count"]
        total_time = agg["total_time"]        # total time in seconds
        avg_time = total_time / count if count > 0 else 0.0
        total_bytes = agg["total_bytes"]
        total_sp = agg["total_sp"]
        total_dp = agg["total_dp"]
        total_hp = agg["total_hp"]

        # Calculate arithmetic intensity (FLOPs/Byte) for each precision
        if total_bytes > 0:
            ai_sp = total_sp / total_bytes
            ai_dp = total_dp / total_bytes
            ai_hp = total_hp / total_bytes
        else:
            # Avoid division by zero: if no bytes, set AI to 0 (or inf if flops exist with zero bytes)
            ai_sp = float('inf') if total_sp > 0 else 0.0
            ai_dp = float('inf') if total_dp > 0 else 0.0
            ai_hp = float('inf') if total_hp > 0 else 0.0

        # Calculate achieved performance (FLOPs/second) for each precision
        if total_time > 0:
            perf_sp = 1e-12 * total_sp / total_time
            perf_dp = 1e-12 * total_dp / total_time
            perf_hp = 1e-12 * total_hp / total_time
        else:
            perf_sp = perf_dp = perf_hp = 0.0

        writer.writerow([
            kname,
            count,
            f"{avg_time:.6f}",
            f"{total_bytes:.0f}",
            total_sp,
            total_dp,
            total_hp,
            f"{ai_sp:.6f}",
            f"{ai_dp:.6f}",
            f"{ai_hp:.6f}",
            f"{perf_sp:.6f}",
            f"{perf_dp:.6f}",
            f"{perf_hp:.6f}"
        ])



# Step 2: Define theoretical hardware ceilings for A100 (SXM4 40GB)&#8203;:contentReference[oaicite:9]{index=9}
peak_mem_bandwidth = 1.555e12  # 1555 GB/s in bytes/sec&#8203;:contentReference[oaicite:10]{index=10}
peak_dp = 9.75e12   # 9.75 TFLOP/s peak FP64&#8203;:contentReference[oaicite:11]{index=11}
peak_sp = 1.95e13   # 19.5 TFLOP/s peak FP32&#8203;:contentReference[oaicite:12]{index=12}
peak_hp = 3.12e14   # 312  TFLOP/s peak FP16 (TensorCore)&#8203;:contentReference[oaicite:13]{index=13}

# Step 1: Load summarized metrics data (AI and FLOPs/s for each kernel)
df = pd.read_csv('roofline_metrics.csv')
  
# First, select entries where at least one AI is greater than zero
df = df[
    (df['AI (DP)'] > 0) | 
    (df['AI (SP)'] > 0) | 
    (df['AI (HP)'] > 0)
]




# Step 3: Initialize Dash app and layout
app = dash.Dash(__name__)
app.layout = html.Div([
    html.H3("GPU Roofline Analysis (A100)", style={'padding-bottom': '10px'}),
    dcc.Checklist(
        id='precision-checklist',
        options=[
            {'label': html.Span('Double-Precision (FP64)', style={'padding-right': '20px'}), 'value': 'DP'},
            {'label': html.Span('Single-Precision (FP32)', style={'padding-right': '20px'}), 'value': 'SP'},
            {'label': html.Span('Half-Precision (FP16)', style={'padding-right': '20px'}), 'value': 'HP'}
        ],
        value=['DP'],
        inline=True,
        style={'padding-bottom': '10px'}
    ),

    html.Div([
        dcc.Graph(id='roofline-graph'),
    ], style={
        'border': '1px solid #ccc',
        'padding': '10px',
        'margin-bottom': '10px',
        'border-radius': '5px'
    }),

    html.Div([
        dcc.Graph(id='dp-ai-histogram'),
        html.Div([
            html.Label('Histogram x_bin_num:', style={'padding-right': '10px'}),
            dcc.Input(id='histogram-x-bin-num', type='number', value=50, min=1, step=1)
        ], style={'padding-top': '10px', 'display': 'flex', 'align-items': 'center'}),
    ], style={
        'border': '1px solid #ccc',
        'padding': '10px',
        'margin-bottom': '10px',
        'border-radius': '5px'
    }),

    html.Div([
        dcc.Graph(id='efficiency-plot'),
        html.Div([
            html.Label('Heatmap x_bin_num:', style={'padding-right': '10px'}),
            dcc.Input(id='heatmap-x-bin-num', type='number', value=10, min=1, step=1),
            html.Label('Heatmap y_bin_num:', style={'padding': '0 10px 0 20px'}),
            dcc.Input(id='heatmap-y-bin-num', type='number', value=10, min=1, step=1)
        ], style={'padding-top': '10px', 'display': 'flex', 'align-items': 'center'}),
    ], style={
        'border': '1px solid #ccc',
        'padding': '10px',
        'margin-bottom': '10px',
        'border-radius': '5px'
    }),
], style={'width': '97%', 'margin': 'auto', 'font-family': 'Arial, sans-serif'})


# Step 4: Define callback to update graph when SP/HP are toggled
@app.callback(
    dash.dependencies.Output('roofline-graph', 'figure'),
    [dash.dependencies.Input('precision-checklist', 'value')]
)
def update_roofline(selected_precisions):

    # Kernel metadata
    kernel_names = df['Kernel Name']
    kernel_times = df['Average Duration (s)']

    # Final arrays
    dp_ai = df['AI (DP)'];    dp_perf = df['Achieved Performance (DP)']
    sp_ai = df['AI (SP)'];    sp_perf = df['Achieved Performance (SP)']
    hp_ai = df['AI (HP)'];    hp_perf = df['Achieved Performance (HP)']

    # Compute the "knee" point (AI where compute = memory limit) for each precision&#8203;:contentReference[oaicite:14]{index=14}
    dp_knee_AI = peak_dp / peak_mem_bandwidth
    sp_knee_AI = peak_sp / peak_mem_bandwidth
    hp_knee_AI = peak_hp / peak_mem_bandwidth

    # Choose a small AI for line start and a max AI for line end (cover data range)
    min_ai = 1e-6#max(1e-3, df[['AI (DP)','AI (SP)','AI (HP)']].min().min())  # smallest non-zero AI in data or ~1e-3
    max_ai = 1e3#max(df[['AI (DP)','AI (SP)','AI (HP)']].max().max(), hp_knee_AI * 10)  # max data AI or a bit beyond HP knee
    # (Using hp_knee*10 ensures we show some flat portion even if data AI is low&#8203;:contentReference[oaicite:15]{index=15})

    # Generate coordinates for each precision's roofline:
    # Memory-bound segment: from min_ai to AI_knee (P = BW * AI)
    # Compute-bound segment: from AI_knee to max_ai (P = peak compute)
    dp_mem_x  = [min_ai, dp_knee_AI];         dp_mem_y  = [1e-12 * min_ai * peak_mem_bandwidth, 1e-12 * peak_dp]
    dp_comp_x = [dp_knee_AI, max_ai];         dp_comp_y = [1e-12 * peak_dp, 1e-12 * peak_dp]
    sp_mem_x  = [min_ai, sp_knee_AI];         sp_mem_y  = [1e-12 * min_ai * peak_mem_bandwidth, 1e-12 * peak_sp]
    sp_comp_x = [sp_knee_AI, max_ai];         sp_comp_y = [1e-12 * peak_sp, 1e-12 * peak_sp]
    hp_mem_x  = [min_ai, hp_knee_AI];         hp_mem_y  = [1e-12 * min_ai * peak_mem_bandwidth, 1e-12 * peak_hp]
    hp_comp_x = [hp_knee_AI, max_ai];         hp_comp_y = [1e-12 * peak_hp, 1e-12 * peak_hp]

    traces = []

    # If HP selected, add HP points and roofline
    if 'DP' in selected_precisions or len(selected_precisions) == 0:
      traces.extend([
          go.Scatter(x=dp_ai, y=dp_perf, mode='markers', name='DP Kernels',
                    marker=dict(color='blue', symbol='circle', size=8),
                    text=kernel_names,
                    hovertemplate=("Kernel: %{text}<br>AI: %{x:.3f} FLOPs/Byte"
                                    "<br>Performance: %{y:.3e} FLOPs/s<extra></extra>")),
          go.Scatter(x=dp_mem_x, y=dp_mem_y, mode='lines', name=f"DP Memory (1555 GB/s)",
                    line=dict(color='blue', width=2, dash='dash'),
                    hovertemplate="DP memory-bound: 1.555e12 * AI<extra></extra>"),
          go.Scatter(x=dp_comp_x, y=dp_comp_y, mode='lines', name=f"DP Compute (Peak 9.75 TFLOP/s)",
                    line=dict(color='blue', width=2, dash='solid'),
                    hovertemplate="DP compute-bound: 9.75e12 FLOPs/s<extra></extra>")
      ])

    # If SP selected, add SP points and roofline
    if 'SP' in selected_precisions:
        # Filter out any zero values (if any) before plotting on log scale
        sp_mask = (sp_ai > 0) & (sp_perf > 0)
        traces.extend([
            go.Scatter(x=sp_ai[sp_mask], y=sp_perf[sp_mask], mode='markers', name='SP Kernels',
                       marker=dict(color='orange', symbol='square', size=7),
                       text=kernel_names[sp_mask],
                       hovertemplate=("Kernel: %{text}<br>AI: %{x:.3f} FLOPs/Byte"
                                      "<br>Performance: %{y:.3e} FLOPs/s<extra></extra>")),
            go.Scatter(x=sp_mem_x, y=sp_mem_y, mode='lines', name=f"SP Memory (1555 GB/s)",
                       line=dict(color='orange', width=2, dash='dash'),
                       hovertemplate="SP memory-bound: 1.555e12 * AI<extra></extra>"),
            go.Scatter(x=sp_comp_x, y=sp_comp_y, mode='lines', name=f"SP Compute (Peak 19.5 TFLOP/s)",
                       line=dict(color='orange', width=2, dash='solid'),
                       hovertemplate="SP compute-bound: 1.95e13 FLOPs/s<extra></extra>")
        ])

    # If HP selected, add HP points and roofline
    if 'HP' in selected_precisions:
        hp_mask = (hp_ai > 0) & (hp_perf > 0)
        traces.extend([
            go.Scatter(x=hp_ai[hp_mask], y=hp_perf[hp_mask], mode='markers', name='HP Kernels',
                       marker=dict(color='green', symbol='triangle-up', size=7),
                       text=kernel_names[hp_mask],
                       hovertemplate=("Kernel: %{text}<br>AI: %{x:.3f} FLOPs/Byte"
                                      "<br>Performance: %{y:.3e} FLOPs/s<extra></extra>")),
            go.Scatter(x=hp_mem_x, y=hp_mem_y, mode='lines', name=f"HP Memory (1555 GB/s)",
                       line=dict(color='green', width=2, dash='dash'),
                       hovertemplate="HP memory-bound: 1.555e12 * AI<extra></extra>"),
            go.Scatter(x=hp_comp_x, y=hp_comp_y, mode='lines', name=f"HP Compute (Peak 312 TFLOP/s)",
                       line=dict(color='green', width=2, dash='solid'),
                       hovertemplate="HP compute-bound: 3.12e14 FLOPs/s<extra></extra>")
        ])

    # Return the updated figure with the selected traces
    return {
        'data': traces,
        'layout': go.Layout(
            xaxis=dict(title="Arithmetic Intensity (FLOPs/Byte)", type="log"),
            yaxis=dict(title="Achieved Performance (FLOPs/s)", type="log"),
            legend=dict(x=0.01, y=0.99),
            margin=dict(l=60, r=20, t=40, b=60),
            hovermode="closest"
        )
    }

@app.callback(
    dash.dependencies.Output('dp-ai-histogram', 'figure'),
    [dash.dependencies.Input('histogram-x-bin-num', 'value'),
     dash.dependencies.Input('precision-checklist', 'value')]
)
def update_histogram(x_bin_num,selected_precisions):

    selected_ai_list = []

    if len(selected_precisions) == 0 or 'DP' in selected_precisions:
        dp_ai_positive = df.loc[
            (df['AI (DP)'] > 0) & (df['Achieved Performance (DP)'] > 0),
            'AI (DP)'
        ]
        selected_ai_list.append(dp_ai_positive)

    if 'SP' in selected_precisions:
        sp_ai_positive = df.loc[
            (df['AI (SP)'] > 0) & (df['Achieved Performance (SP)'] > 0),
            'AI (SP)'
        ]
        selected_ai_list.append(sp_ai_positive)

    if 'HP' in selected_precisions:
        hp_ai_positive = df.loc[
            (df['AI (HP)'] > 0) & (df['Achieved Performance (HP)'] > 0),
            'AI (HP)'
        ]
        selected_ai_list.append(hp_ai_positive)

    # Concatenate explicitly positive-only AI values
    ai = pd.concat(selected_ai_list)


    if x_bin_num is None or x_bin_num < 1:
        x_bin_num = 50

    log_bins = np.logspace(np.log10(ai.min()/2.0), np.log10(ai.max()*2.0), num=int(x_bin_num))
    ai_hist, bin_edges = np.histogram(ai, bins=log_bins)
    bin_centers = 10 ** ((np.log10(bin_edges[1:]) + np.log10(bin_edges[:-1])) / 2)

    figure={
        'data': [
            go.Bar(
                x=bin_centers,
                y=ai_hist,
                width=np.diff(bin_edges),
                marker=dict(color='blue'),
                opacity=0.75,
                name="DP AI Histogram"
            )
        ],
        'layout': go.Layout(
            title="Distribution of Arithmetic Intensity (DP Kernels)",
            xaxis=dict(title="Arithmetic Intensity (FLOPs/Byte)", type='log', range=[-6,3]),
            yaxis=dict(title="Kernel Count"),
            margin=dict(l=60, r=20, t=40, b=60)
        )
    }
    return figure

@app.callback(
    dash.dependencies.Output('efficiency-plot', 'figure'),
    [dash.dependencies.Input('heatmap-x-bin-num', 'value'),
     dash.dependencies.Input('heatmap-y-bin-num', 'value'),
     dash.dependencies.Input('precision-checklist', 'value')]
)
def update_efficiency_plot(x_bin_num, y_bin_num, selected_precisions):

    selected_ai_list = []
    selected_perf_list = []
    selected_efficiency_list = []
    selected_kernel_times_list = []
    selected_kernel_names_list = []

    # DP precision
    if 'DP' in selected_precisions or len(selected_precisions) == 0:
        dp_valid = df[(df['AI (DP)'] > 0) & (df['Achieved Performance (DP)'] > 0)]
        dp_ai_values = dp_valid['AI (DP)']
        dp_perf_values = dp_valid['Achieved Performance (DP)']
        dp_kernel_times = dp_valid['Average Duration (s)']
        dp_kernel_names = dp_valid['Kernel Name']
        theoretical_max_dp = np.minimum(1e-12 * peak_dp, 1e-12 * dp_ai_values * peak_mem_bandwidth)
        efficiency_dp = dp_perf_values / theoretical_max_dp

        selected_ai_list.append(dp_ai_values)
        selected_perf_list.append(dp_perf_values)
        selected_efficiency_list.append(efficiency_dp)
        selected_kernel_times_list.append(dp_kernel_times)
        selected_kernel_names_list.append(dp_kernel_names)

    # SP precision
    if 'SP' in selected_precisions:
        sp_valid = df[(df['AI (SP)'] > 0) & (df['Achieved Performance (SP)'] > 0)]
        sp_ai_values = sp_valid['AI (SP)']
        sp_perf_values = sp_valid['Achieved Performance (SP)']
        sp_kernel_times = sp_valid['Average Duration (s)']
        sp_kernel_names = sp_valid['Kernel Name']
        theoretical_max_sp = np.minimum(1e-12 * peak_sp, 1e-12 * sp_ai_values * peak_mem_bandwidth)
        efficiency_sp = sp_perf_values / theoretical_max_sp

        selected_ai_list.append(sp_ai_values)
        selected_perf_list.append(sp_perf_values)
        selected_efficiency_list.append(efficiency_sp)
        selected_kernel_times_list.append(sp_kernel_times)
        selected_kernel_names_list.append(sp_kernel_names)

    # HP precision
    if 'HP' in selected_precisions:
        hp_valid = df[(df['AI (HP)'] > 0) & (df['Achieved Performance (HP)'] > 0)]
        hp_ai_values = hp_valid['AI (HP)']
        hp_perf_values = hp_valid['Achieved Performance (HP)']
        hp_kernel_times = hp_valid['Average Duration (s)']
        hp_kernel_names = hp_valid['Kernel Name']
        theoretical_max_hp = np.minimum(1e-12 * peak_hp, 1e-12 * hp_ai_values * peak_mem_bandwidth)
        efficiency_hp = hp_perf_values / theoretical_max_hp

        selected_ai_list.append(hp_ai_values)
        selected_perf_list.append(hp_perf_values)
        selected_efficiency_list.append(efficiency_hp)
        selected_kernel_times_list.append(hp_kernel_times)
        selected_kernel_names_list.append(hp_kernel_names)

    # Concatenate all collected arrays
    ai = pd.concat(selected_ai_list)
    perf = pd.concat(selected_perf_list)
    efficiency = pd.concat(selected_efficiency_list)
    kernel_times = pd.concat(selected_kernel_times_list)
    kernel_names = pd.concat(selected_kernel_names_list)


    if x_bin_num is None or x_bin_num < 1:
        x_bin_num = 10

    if y_bin_num is None or y_bin_num < 1:
        y_bin_num = 10

    #theoretical_max = np.minimum(1e-12 * peak_dp, 1e-12 * ai * peak_mem_bandwidth)
    #efficiency = perf / theoretical_max

    efficiency = np.clip(efficiency, 0, 1)

    x_bins = np.logspace(np.log10(ai.min()/2), np.log10(ai.max()*2), num=int(x_bin_num)+1)
    y_bins = np.linspace(0, 1.0, int(y_bin_num)+1)

    heatmap_data = np.zeros((len(y_bins)-1, len(x_bins)-1))

    #print(ai['Kernel Name'])

    for i in range(len(ai)):

        x_idx = np.digitize(ai.iloc[i], x_bins, right=True) - 1
        y_idx = np.digitize(efficiency.iloc[i], y_bins, right=True) - 1

        if 0 <= x_idx < heatmap_data.shape[1] and 0 <= y_idx < heatmap_data.shape[0]:
            heatmap_data[y_idx, x_idx] += kernel_times.values[i]

    if heatmap_data.sum() > 0:
        heatmap_data /= heatmap_data.sum()

    figure={
        'data': [
            go.Heatmap(
                z=heatmap_data,
                x=x_bins,
                y=y_bins,
                colorscale=[
                    [0.0, 'rgba(255,255,255,0)'],
                    [1e-9, 'rgb(255,245,240)'],
                    [1.0, 'rgb(150,0,13)']
                ],
                colorbar=dict(
                    x=0.95,
                    xanchor='right',
                    y=0.5,
                    yanchor='middle',
                    len=0.6
                ),
                zmin=0, zmax=heatmap_data.max(),
                showscale=True,
                hovertemplate=(
                    "AI bin: %{x:.2f}<br>Efficiency bin: %{y:.2f}<br>"
                    "Time fraction: %{z:.2%}<extra></extra>"
                )
            ),
            go.Scatter(
                x=ai,
                y=efficiency,
                mode='markers',
                name='Efficiency',
                marker=dict(color='green', size=8),
                text=kernel_names,
                hovertemplate=("Kernel: %{text}<br>AI: %{x:.3f} FLOPs/Byte"
                               "<br>Efficiency: %{y:.2f}<extra></extra>")
            )
        ],
        'layout': go.Layout(
            title="DP Kernel Efficiency Relative to Roofline",
            xaxis=dict(title="Arithmetic Intensity (FLOPs/Byte)", type="log", range=[-6,3]),
            yaxis=dict(title="Efficiency (Achieved / Theoretical)", range=[0, 1.05]),
            margin=dict(l=60, r=20, t=40, b=60),
            hovermode="closest"
        )
    }
    return figure


# Step 5: Launch the Dash app server (opens the dashboard in a web browser)
if __name__ == '__main__':
    app.run_server(debug=True)
