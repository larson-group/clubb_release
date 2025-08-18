
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
totals_df['Vectorization Ratio'] = totals_df['FP_ARITH_INST_RETIRED.VECTOR'] / (
    totals_df['FP_ARITH_INST_RETIRED.VECTOR'] + totals_df['FP_ARITH_INST_RETIRED.SCALAR']
)


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
                            name='Vectorization Ratio'
                        )
                    ],
                    'layout': go.Layout(
                        title=dict(
                            text="CLUBB Vectorization Ratio on Intel 6430 CPU",
                            font=dict(size=24)  # Title font size
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
            )
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
                options=[{'label': col, 'value': col} for col in ['CPI'] + [
                    col for col in totals_df.columns if col not in ['CPI', 'Vectorization Ratio']]],
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
    for metric in selected_metrics:
        #y_safe = [max(y, 1) for y in totals_df[metric]]
        y_safe = totals_df[metric]
        fig.add_trace(go.Scatter(
            x=totals_df.index,
            y=y_safe,  # Use the cleaned data
            mode='lines+markers',
            name=metric
        ))
    fig.update_layout(
        title=dict(
            text="CLUBB CPI on Intel 6430 CPU",
            font=dict(size=24)  # Title font size
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
            type="log",
            #dtick=4
        ),
        yaxis=dict(
            title=dict(
                text="Cycles per Instruction (CPI)",
                font=dict(size=18)
            ),
            showline=True,
            linecolor="black",
            linewidth=1,
            mirror=True,
            #type="log",
            #dtick=0.1
        ),
        margin=dict(l=10, r=10, t=50, b=10),
        template="plotly_white"
    )
    return fig

if __name__ == '__main__':
    app.run(debug=True)
