#!/bin/python3

import os
import pandas as pd
import re
from glob import glob
from dash import Dash, dcc, html, Input, Output, State, dash_table
import plotly.express as px
import numpy as np
from collections import defaultdict
from dash.dependencies import Input, Output, ALL
from dash.exceptions import PreventUpdate
import pyperclip
from scipy.optimize import minimize

# Functions safe to use in the custom plot
safe_functions = {
    "abs": abs,
    "max": max,
    "min": min,
    "round": round,
    "ceil": np.ceil,
    "floor": np.floor,
    "sqrt": np.sqrt,
    "log": np.log,
    "log2": np.log2,
    "log10": np.log10,
    "exp": np.exp,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "pi": np.pi,
    "e": np.e,
    "np": np
}

gpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m", "id": "m_val"},
    {"name": "b", "id": "b_val"},
    {"name": "b_est", "id": "b_est_val"},
    {"name": "RMSE (abs '%' diff)", "id": "rms_error"},
]

cpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m", "id": "m_val"},
    {"name": "m_est", "id": "m_est_val"},
    {"name": "b", "id": "b_val"},
    {"name": "b_est", "id": "b_est_val"},
    {"name": "c", "id": "c_val"},
    {"name": "k", "id": "k_val"},
    {"name": "o", "id": "o_val"},
    {"name": "Cache Pen Func", "id": "cp_func"},
    {"name": "RMSE (abs '%' diff)", "id": "rms_error"},
]

vcpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "T_v", "id": "T_v_val"},
    {"name": "T_r", "id": "T_r_val"},
    {"name": "b", "id": "b_val"},
    {"name": "b_est", "id": "b_est_val"},
    {"name": "c", "id": "c_val"},
    {"name": "k", "id": "k_val"},
    {"name": "o", "id": "o_val"},
    {"name": "Cache Pen Func", "id": "cp_func"},
    {"name": "RMSE (abs '%' diff)", "id": "rms_error"},
]


plot_div_style = {
    "margin-bottom": "2%",
    "border": "1px solid black",
    "padding": "1%",
    "width": "auto",              
    "maxWidth": "100%",           # Prevent overflow
    "height": "35%",             # Let it fill available vertical space
    "flexGrow": 1,                # Tell flex layout to expand this div if space allows
    "minWidth": 0,                # Prevent layout overflow in flexbox
    "display": "flex",            # Optional, useful if you're aligning internal content
    "flexDirection": "column",    # Keeps label/table below the graph
    "alignItems": "center"
}

fit_plot_div_style = {
    "margin-bottom": "2%",
    "border": "1px solid black",
    "padding": "1%",
    "width": "auto",              
    "maxWidth": "100%",           # Prevent overflow
    "height": "60%",             # Let it fill available vertical space
    "flexGrow": 1,                # Tell flex layout to expand this div if space allows
    "minWidth": 0,                # Prevent layout overflow in flexbox
    "display": "flex",            # Optional, useful if you're aligning internal content
    "flexDirection": "column",    # Keeps label/table below the graph
    "alignItems": "left"
}

graph_style = {"width": "100%", "height": "100%"}

cache_pen_funcs = [ "loglog" ] #, "sigmoid", "sqrtlog" ]

def model_vcpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    cols_per_core = ngrdcol / N_tasks
    flops_per_vop = N_vsize / N_prec

    b, c, k, o = params

    T_r = runtime[np.where(cols_per_core == 1)] - b
    T_v = runtime[np.where(cols_per_core == flops_per_vop)] - b

    N_s = np.ceil( cols_per_core ) % ( flops_per_vop )
    N_v = np.floor( cols_per_core / ( flops_per_vop ) )

    f = b + T_r * N_s + T_v * N_v

    problem_size = (ngrdcol * N_vlevs / 1000.0) * (N_prec / 64)

    if cp_func == "sigmoid":
        p = 1 + c / (1 + np.exp(-k * (problem_size - o)))
    elif cp_func == "loglog":
        p =  c * np.log(np.log( np.exp((problem_size - o)) + 1 ) + 1 ) + 1
    elif cp_func == "sqrtlog":
        p = c * np.sqrt( np.log( np.exp(k*(problem_size - o)) + 1 ) ) + 1
    else:
        p = 1

    T_cpu = p * f

    return T_cpu


def model_cpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    m, b, c, k, o = params

    f = b + m * ngrdcol

    problem_size = (ngrdcol * N_vlevs / 1000.0) * (N_prec / 64)

    if cp_func == "sigmoid":
        p = 1 + c / (1 + np.exp(-k * (problem_size - o )))
    elif cp_func == "loglog":
        p =  c * np.log(np.log( np.exp((problem_size - o)) + 1 ) + 1 ) + 1
    elif cp_func == "sqrtlog":
        p = c * np.sqrt( np.log( np.exp(k*(problem_size - o)) + 1 ) ) + 1
    else:
        p = 1

    T_cpu = p * f

    return T_cpu


def model_gpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec):

    m, b = params

    f = b + m * ngrdcol

    T_gpu = f

    return T_gpu


def vcpu_objective(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):

    T_cpu = model_vcpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )

    # rms_error = np.sqrt( np.mean( ( ngrdcol / runtime - ngrdcol / T_cpu )**2 ) ) \
    #             / ( np.max(ngrdcol / runtime) - np.min(ngrdcol / runtime) )

    cps = ngrdcol / runtime
    cps_model = ngrdcol / T_cpu

    rms_error = np.sqrt( np.mean( ( ( cps - cps_model ) / cps * 100 )**2 ) )

    return rms_error

    
def cpu_objective(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):

    T_cpu = model_cpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )

    # rms_error = np.sqrt( np.mean( ( ngrdcol / runtime - ngrdcol / T_cpu )**2 ) ) \
    #             / ( np.max(ngrdcol / runtime) - np.min(ngrdcol / runtime) )

    cols_per_core = ngrdcol / N_tasks
    flops_per_vop = N_vsize / N_prec

    cps = ngrdcol / runtime
    cps_model = ngrdcol / T_cpu

    abs_percent_diff = ( cps - cps_model ) / cps

    rms_error = np.sqrt( np.mean( ( abs_percent_diff[np.where(cols_per_core % flops_per_vop == 0)] * 100 )**2 ) )

    return rms_error


def gpu_objective(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec ):

    T_gpu = model_gpu_time(ngrdcol, runtime, params, N_tasks, N_vsize, N_prec )

    cps = ngrdcol / runtime
    cps_model = ngrdcol / T_gpu

    #abs_percent_diff = ( cps - cps_model ) / cps
    #rms_error = np.sqrt( np.mean( ( abs_percent_diff[np.where(ngrdcol >= 10000)] * 100 )**2 ) )

    abs_diff = np.abs( cps - cps_model )
    rms_error = np.sqrt( np.mean( abs_diff**2 ) )

    return rms_error


def model_throughputs(df, column_name, N_tasks, N_vsize, N_prec, N_vlevs, model_version):

    ngrdcol = df["ngrdcol"].values  # Access ngrdcol if needed
    runtime = df[column_name].values  # Convert selected column to numpy array


    #============================== VCPU Model ==============================

    if model_version == "vcpu":

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]

        initial_guess = [b_est, 3.0, 0.02, 300]

        cache_pen_best = None
        rms_min = None
        params_opt = None

        for cp_func in cache_pen_funcs:

            result = minimize(  vcpu_objective, initial_guess, 
                                args=(ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func), 
                                method='Nelder-Mead')

            rms_error = result.fun

            if rms_min is None or rms_error < rms_min:
                cache_pen_best = cp_func
                rms_min = rms_error
                params_opt = result.x


        T_cpu = model_vcpu_time(ngrdcol, runtime, params_opt, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        b = params_opt[0]
        T_r = runtime[np.where(cols_per_core == 1)] - b
        T_v = runtime[np.where(cols_per_core == flops_per_vop)] - b

        fit_params = {
            "c_val": params_opt[1],
            "k_val": params_opt[2],
            "o_val": params_opt[3], 
            "b_val": b, 
            "b_est_val": b_est, 
            "T_v_val": T_v[0], 
            "T_r_val": T_r[0], 
            "cp_func": cache_pen_best, 
            "rms_error" : rms_error
        }

    #============================== CPU Model ==============================

    if model_version == "cpu":

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]
        m_est = ( runtime[3] - b_est ) / ( ngrdcol[3] )

        initial_guess = [m_est, b_est, 3.0, 0.02, 300]

        cache_pen_best = None
        rms_min = None
        params_opt = None

        for cp_func in cache_pen_funcs:

            result = minimize(  cpu_objective, initial_guess, 
                                args=(ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func), 
                                method='Nelder-Mead')

            rms_error = result.fun

            if rms_min is None or rms_error < rms_min:
                cache_pen_best = cp_func
                rms_min = rms_error
                params_opt = result.x


        T_cpu = model_cpu_time(ngrdcol, runtime, params_opt, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        #params_opt = np.exp( params_opt_scaled )

        fit_params = {
            "m_val": params_opt[0], 
            "m_est_val": m_est, 
            "b_val": params_opt[1], 
            "b_est_val": b_est, 
            "c_val": params_opt[2],
            "k_val": params_opt[3],
            "o_val": params_opt[4], 
            "cp_func": cache_pen_best, 
            "rms_error" : rms_error
        }

    #============================== GPU Model ==============================

    if model_version == "gpu":

        b_est = runtime[-2] - ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] ) * ngrdcol[-2]
        m_est = ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] )

        initial_guess = [.1, .01]

        rms_min = None
        params_opt = None


        result = minimize(  gpu_objective, initial_guess, 
                            args=(ngrdcol, runtime, N_tasks, N_vsize, N_prec), 
                            method='Nelder-Mead')

        rms_error = result.fun
        params_opt = result.x

        T_cpu = model_gpu_time(ngrdcol, runtime, params_opt, N_tasks, N_vsize, N_prec)

        fit_params = {
            "m_val": params_opt[0], 
            "m_est_val": m_est, 
            "b_val": params_opt[1], 
            "b_est_val": b_est, 
            "rms_error" : rms_error
        }


    return pd.Series(T_cpu, index=df.index), fit_params


def get_shared_variables(csv_files):
    variable_sets = {}
    for file in csv_files:
        try:
            with open(file, 'r') as f:
                header = f.readline().strip()
                variables = header.split(",")[1:]
                variable_sets[file] = set(variables)
        except Exception as e:
            print(f"Error reading {file}: {e}")

    if not variable_sets:
        return [], variable_sets

    shared_vars = sorted(set.intersection(*variable_sets.values())) if variable_sets.values() else []
    return shared_vars, variable_sets

def get_all_variables(csv_files):
    variable_sets = {}
    all_variables = set()
    
    for file in csv_files:
        try:
            with open(file, 'r') as f:
                header = f.readline().strip()
                variables = header.split(",")[1:]
                variable_sets[file] = set(variables)
                all_variables.update(variables)
        except Exception as e:
            print(f"Error reading {file}: {e}")

    return sorted(all_variables), variable_sets

def natural_key(string):
    """Helper function to split a string into a list of numbers and text for natural sorting."""
    return [int(part) if part.isdigit() else part for part in re.split(r'(\d+)', string)]

def group_files_by_case(csv_files):
    cases = defaultdict(dict)
    for file in csv_files:
        match = re.match(r"^(.*)_(\w+)\.csv$", os.path.basename(file))
        if match:
            filename, case = match.groups()
            cases[case][filename] = file

    # Sort the filenames within each case
    sorted_cases = {}
    for case, files in cases.items():
        sorted_cases[case] = dict(sorted(files.items(), key=lambda item: natural_key(item[0])))

    return sorted_cases

def evaluate_function(func, df):
    try:
        local_vars = {col: df[col] for col in df.columns}
        return eval(func, {"__builtins__": None}, {**local_vars, **safe_functions})
    except Exception as e:
        print(f"Error evaluating function '{func}': {e}")
        return None

def plot_with_enhancements(fig, title):
    fig.update_traces(mode="lines+markers")
    fig.update_layout(
        title=title,
        xaxis=dict(showgrid=True, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridcolor="lightgray", zeroline=True, range=[0, None]),
        margin=dict(l=10, r=10, t=30, b=10),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    fig.update_xaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    return fig

def launch_dash_app(grouped_files, all_variables):

    data = {case: {filename: pd.read_csv(filepath, comment="#") for filename, filepath in files.items()} for case, files in grouped_files.items()}

    app = Dash(__name__)
    app.title = "Dynamic Plotter"
                   
    app.layout = html.Div([
        html.Div([
            html.Div(dcc.Graph(id="plot-columns-per-second", config={"responsive": True}, style=graph_style), style=plot_div_style),
            html.Div(dcc.Graph(id="plot-raw", config={"responsive": True}, style=graph_style), style=plot_div_style),
            html.Div([
                dcc.Graph(id="fit-plot", config={"responsive": True}, style=graph_style), 
                html.Label("Vector CPU Model"),
                dash_table.DataTable(
                    id='vcpu-param-table',
                    columns=vcpu_model_columns,
                    data=[],  # Initial empty table
                    editable=True,
                    row_deletable=False  # Allows deleting rows directly
                ),
                html.Label("Simple CPU Model"),
                dash_table.DataTable(
                    id='cpu-param-table',
                    columns=cpu_model_columns,
                    data=[],  # Initial empty table
                    editable=True,
                    row_deletable=False  # Allows deleting rows directly
                ),
                html.Label("GPU Model"),
                dash_table.DataTable(
                    id='gpu-param-table',
                    columns=gpu_model_columns,
                    data=[],  # Initial empty table
                    editable=True,
                    row_deletable=False  # Allows deleting rows directly
                ),
            ], style=fit_plot_div_style),
            html.Div([
                dcc.Graph(id="plot-custom", config={"responsive": True}, style=graph_style),
                dcc.Input(
                    id="custom-function",
                    type="text",
                    placeholder="Enter function (e.g., baseline + 5 * iter_total)",
                    debounce=True,
                    style={"width": "30%", "margin-top": "10px", "text-align": "center"}
                ),
            ], style=plot_div_style),
        ], style={  "flex": "1", 
                    "padding-right": "10px", 
                    "align-items": "center",
                    "width": "auto",
                    "height": "auto" }
        ),

        html.Div([
            html.H4("Select CSV files to display:"),
            html.Div([
                html.Details(
                    [
                        html.Summary(case),
                        dcc.Checklist(
                            id={"type": "file-checkbox", "case": case},
                            options=[
                                {"label": filename, "value": f"{case}/{filename}"}
                                for filename in files.keys()
                            ],
                            value=[
                                f"{case}/{filename}"
                                for filename in list(files.keys())[:2]
                            ] if i == 0 and files else []
                        )
                    ],
                    open=(i == 0)
                )
                for i, (case, files) in enumerate(grouped_files.items())
            ], style={"margin-bottom": "20px"}),

            html.H4("Select variable to plot:"),
            dcc.Dropdown(
                id="variable-dropdown",
                options=[{"label": var, "value": var} for var in all_variables],
                value="compute_i",
                style={"margin-bottom": "20px"}
            ),

            html.Button(
                "Copy CSV Data to Clipboard",
                id="copy-csv-button",
                n_clicks=0,
                style={"margin-top": "10px", "padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
            ),

            html.Div([
                html.Label("x-axis:", style={"margin-right": "10px"}),
                dcc.RadioItems(
                    id="x-axis-scale",
                    options=[
                        {"label": "Linear", "value": "linear"},
                        {"label": "Log", "value": "log"}
                    ],
                    value="log",
                    inline=True,
                )
            ], style={"display": "flex", "align-items": "center", "margin-top": "10px"}),

            html.Div([
                html.Label("y-axis:", style={"margin-right": "10px"}),
                dcc.RadioItems(
                    id="y-axis-scale",
                    options=[
                        {"label": "Linear", "value": "linear"},
                        {"label": "Log", "value": "log"}
                    ],
                    value="linear",
                    inline=True,
                )
            ], style={"display": "flex", "align-items": "center", "margin-top": "10px"}),

            # New div for the fit button and input fields
            html.Div([
                html.Label("Fit Function Variables:", style={"font-weight": "bold", "margin-bottom": "5px"}),

                # Stack inputs vertically
                html.Div([
                    dcc.RadioItems(
                        id="fit-plot-mode",
                        options=[
                            {"label": "Time", "value": "time"},
                            {"label": "Cols per Sec", "value": "cps"}
                        ],
                        value="cps",
                        inline=True,
                    ),
                    #html.Label("Total Cores:"),
                    #dcc.Input(id="N_tasks", type="number", value=128, style={"width": "60%"}),
                    html.Label("Vector Length:", style={"margin-top": "5px"}),
                    dcc.Input(id="N_vsize", type="number", value=256, style={"width": "60%"}),
                    html.Label("Floating Point Precision:", style={"margin-top": "5px"}),
                    dcc.Input(id="N_prec", type="number", value=64, style={"width": "60%"}),
                ], style={"display": "flex", "flex-direction": "column", "gap": "5px", "margin-bottom": "10px"}),

                html.Button(
                    "Fit",
                    id="fit-function-button",
                    n_clicks=0,
                    style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
                ),
            ], style={"margin-top": "10px", "padding": "10px", "border": "1px solid #000", "background-color": "#f9f9f9"}),

        ], style={
            "width": "20%",
            "padding": "10px",
            "border": "1px solid black",
            "position": "sticky",
            "top": "0",
            "height": "fit-content",
            "background-color": "white",
            "z-index": "1000"
        })
    ], style={"display": "flex", "flexDirection": "row", "width": "100%", "height": "100vh"})

    @app.callback(
    Output("plot-raw", "figure"),
        [
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value")
        ]
    )
    def update_raw_plot(selected_files, selected_variable, xaxis_scale, yaxis_scale):
        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()
        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:
                df = data[case][filename]
                if selected_variable in df.columns:  # Check before accessing
                    temp_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                    temp_df["Source"] = f"{case}/{filename}"
                    combined_df = pd.concat([combined_df, temp_df])
        fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Source")
        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            autosize=True,
            uirevision='constant'
        )
        return plot_with_enhancements(fig, f"Runtime of '{selected_variable}' vs. Number of Grid Columns")

    @app.callback(
    Output("plot-columns-per-second", "figure"),
        [
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value")  
        ]
    )
    def update_columns_per_second_plot(selected_files, selected_variable, xaxis_scale, yaxis_scale):

        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()
        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:
                df = data[case][filename]
                if selected_variable in df.columns:  # Check before accessing
                    temp_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                    temp_df["Source"] = f"{case}/{filename}"
                    temp_df["Columns per Second"] = temp_df["ngrdcol"] / temp_df[selected_variable]
                    combined_df = pd.concat([combined_df, temp_df])
        fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Source")
        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            autosize=True,
            uirevision='constant'
        )
        return plot_with_enhancements(fig, f"Throughput (ngrdcol/{selected_variable}) vs. Number of Grid Columns")

    @app.callback(
    Output("plot-custom", "figure"),
        [
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("custom-function", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value")  
        ]
    )
    def update_custom_plot(selected_files, custom_function, xaxis_scale, yaxis_scale):

        if not selected_files or not custom_function:
            return px.scatter(title="Enter a custom function and select CSV files")

        # Flatten the list of selected files
        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()

        for case, filename in flat_files:
            if case in data and filename in data[case]:
                temp_df = data[case][filename].copy()
                temp_df["Source"] = f"{case}/{filename}"

                # Evaluate the custom function
                result = evaluate_function(custom_function, temp_df)
                if result is not None:
                    temp_df["Custom Function"] = result
                    combined_df = pd.concat([combined_df, temp_df])

        if combined_df.empty:
            return px.scatter(title="No valid data for the custom function")

        fig = px.line(combined_df, x="ngrdcol", y="Custom Function", color="Source")
        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            autosize=True,
            uirevision='constant'
        )
        return plot_with_enhancements(fig, "Custom Function vs. Number of Grid Columns")

    @app.callback(
    Output("copy-csv-button", "children"),  # Update button text to confirm copy
        [
            Input("copy-csv-button", "n_clicks"),
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value")
        ]
    )
    def copy_csv_to_clipboard(n_clicks, selected_files, selected_variable):

        if not n_clicks or not selected_files:
            raise PreventUpdate

        # Flatten the list of selected files: list of (case, filename) pairs
        flat_files = [item.split("/") for sublist in selected_files for item in sublist]

        # The order in 'header' will match the order of these files
        header = ["ngrdcol"]
        header += [f"{case}/{filename}" for (case, filename) in flat_files]

        # This will map: ngrdcol -> { "case1/file1": val, "case2/file2": val, ... }
        data_rows = defaultdict(dict)

        # Populate data_rows
        for (case, filename) in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:
                df = data[case][filename][["ngrdcol", selected_variable]].copy()
                for ngrdcol_val in df["ngrdcol"].unique():
                    if ngrdcol_val not in data_rows:
                        data_rows[ngrdcol_val] = {}
                    for idx, row in df[df["ngrdcol"] == ngrdcol_val].iterrows():
                        data_rows[ngrdcol_val][f"{case}/{filename}"] = row.get(selected_variable, "")

        # Rebuild and sort the union of all ngrdcol values
        all_ngrdcols = sorted(set(data_rows.keys()))

        # Build CSV lines
        csv_data = []
        csv_data.append(",".join(header))

        for ngrd_val in all_ngrdcols:
            row_values = [str(ngrd_val)]
            # For each file column (in the same order as 'header[1:]'), append its value or blank
            for file_col in header[1:]:
                value = data_rows[ngrd_val].get(file_col, "")  # fill as "" or "NaN"
                row_values.append(str(value))
            csv_data.append(",".join(row_values))

        # Copy CSV data to clipboard
        csv_string = "\n".join(csv_data)
        try:
            pyperclip.copy(csv_string)
            return "Copied to Clipboard!"
        except Exception as e:
            print(f"Error copying to clipboard: {e}")
            return "Copy Failed"

    @app.callback(
        [
            Output("fit-plot", "figure"),
            Output('vcpu-param-table', 'data'),
            Output('cpu-param-table', 'data'),
            Output('gpu-param-table', 'data'),
        ],
        [
            Input("fit-function-button", "n_clicks"),
            #Input("N_tasks", "value"),
            Input("N_vsize", "value"),
            Input("N_prec", "value"),
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("fit-plot-mode", "value")
        ],
    )
    def update_fit_plot(n_clicks, N_vsize, N_prec, selected_files, selected_variable, xaxis_scale, yaxis_scale, plot_mode):
        
        if n_clicks == 0 or not selected_files:
            raise PreventUpdate

        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()

        vcpu_table_data = []
        cpu_table_data = []
        gpu_table_data = []

        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:

                original_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                original_df["Name"] = f"{case}/{filename}"
                original_df["Columns per Second"] = original_df["ngrdcol"] / original_df[selected_variable]

                combined_df = pd.concat([combined_df, original_df])
                
                # Get the number of cores used 
                N_tasks = data[case][filename]["tasks"][0]
                N_vlevs = data[case][filename]["nz"][0]

                # Apply model_throughput to the selected variable and add it as "_dup"
                #temp_df_dup = temp_df.copy()
                vcpu_df = original_df.copy()
                cpu_df = original_df.copy()
                gpu_df = original_df.copy()
                
                if any(gpu_name in filename for gpu_name in ["A100", "V100", "H100"]):

                    gpu_df[selected_variable], gpu_params = model_throughputs(gpu_df, selected_variable, N_tasks, N_vsize, N_prec, N_vlevs, "gpu")
                    gpu_df["Columns per Second"] = gpu_df["ngrdcol"] / gpu_df[selected_variable]
                    gpu_df["Name"] = f"{case}/{filename}_gmodel"

                    gpu_table_data.append({
                        col["id"]: f"{gpu_params.get(col['id']):.3e}" if isinstance(gpu_params.get(col["id"]), (int, float)) else gpu_params.get(col["id"])
                        for col in cpu_model_columns
                    })
                    gpu_table_data[-1]["name"] = f"{case}/{filename}"

                    combined_df = pd.concat([combined_df, gpu_df])

                else:

                    vcpu_df[selected_variable], vcpu_params = model_throughputs(vcpu_df, selected_variable, N_tasks, N_vsize, N_prec, N_vlevs, "vcpu")
                    vcpu_df["Columns per Second"] = vcpu_df["ngrdcol"] / vcpu_df[selected_variable]
                    vcpu_df["Name"] = f"{case}/{filename}_vmodel"

                    vcpu_table_data.append({
                        col["id"]: f"{vcpu_params.get(col['id']):.3e}" if isinstance(vcpu_params.get(col["id"]), (int, float)) else vcpu_params.get(col["id"])
                        for col in vcpu_model_columns
                    })
                    vcpu_table_data[-1]["name"] = f"{case}/{filename}"


                    cpu_df[selected_variable], cpu_params = model_throughputs(cpu_df, selected_variable, N_tasks, N_vsize, N_prec, N_vlevs, "cpu")
                    cpu_df["Columns per Second"] = cpu_df["ngrdcol"] / cpu_df[selected_variable]
                    cpu_df["Name"] = f"{case}/{filename}_cmodel"

                    cpu_table_data.append({
                        col["id"]: f"{cpu_params.get(col['id']):.3e}" if isinstance(cpu_params.get(col["id"]), (int, float)) else cpu_params.get(col["id"])
                        for col in cpu_model_columns
                    })
                    cpu_table_data[-1]["name"] = f"{case}/{filename}"

                    combined_df = pd.concat([combined_df, vcpu_df, cpu_df])
                
        if plot_mode == "cps":
            fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Name")
        else:
            fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Name")

        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            autosize=True,
            uirevision='constant'
        )
        fig = plot_with_enhancements(fig, f"{selected_variable} vs. Number of Grid Columns (Raw Values)")

        return fig.to_dict(), vcpu_table_data, cpu_table_data, gpu_table_data

    app.run(debug=True,port=8051)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot shared variables across multiple CSV files.")
    parser.add_argument("-dir", required=True, help="Directory to look for .csv files.")
    args = parser.parse_args()

    csv_files = glob(os.path.join(args.dir, "*.csv"))
    grouped_files = group_files_by_case(csv_files)

    # shared_variables, variable_sets = get_shared_variables(
    #     [file for files in grouped_files.values() for file in files.values()]
    # )
    all_variables, variable_sets = get_all_variables(
        [file for files in grouped_files.values() for file in files.values()]
    )

    # if not shared_variables:
    #     print("No shared variables found across the CSV files.")
    #     print("Variable presence per file:")

    #     # Get all unique variables across files
    #     all_variables = set.union(*variable_sets.values()) if variable_sets else set()

    #     # Report missing variables for each file
    #     for file, variables in variable_sets.items():
    #         missing_vars = all_variables - variables
    #         if missing_vars:
    #             print(f"{file} is missing the variables: {', '.join(sorted(missing_vars))}")

    #     exit(1)

    #print(grouped_files)

    launch_dash_app(grouped_files, all_variables)
