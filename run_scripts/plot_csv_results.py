#!/bin/python3

import os
import pandas as pd
import re
from glob import glob
from dash import Dash, dcc, html, Input, Output
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

def model_cpu_time(ngrdcol, runtime, params ):

    L = len(runtime)

    N_cores = 128
    N_vsize = 256
    N_prec  = 32

    d = 0
    b = runtime[d] - ( runtime[d+1] - runtime[d] ) / ( ngrdcol[d+1] - ngrdcol[d] ) * ngrdcol[d]

    T_s = runtime[np.where(ngrdcol == N_cores)] - b
    T_v = runtime[np.where(ngrdcol == N_vsize*N_cores/N_prec)] - b

    N_s = np.ceil( ngrdcol / N_cores ) % ( N_vsize / N_prec )
    N_v = np.floor( ngrdcol / ( N_cores * N_vsize / N_prec ) )

    f = b + T_s * N_s + T_v * N_v

    o, k, c = params
    #p = 1 + c / (1 + np.exp(-k * (ngrdcol - o)))
    p = c * np.log( np.log( np.exp(k*(ngrdcol-o)) + 1 ) + 1 ) + 1

    #print(p)
    #print(ngrdcol)

    T_cpu = p * f

    return T_cpu


def objective(params, ngrdcol, runtime):

    T_cpu = model_cpu_time(ngrdcol, runtime, params)

    #N = np.where( ngrdcol == 4096 )[0][0]
    #rms_error = np.sqrt( np.mean( ngrdcol[0:N] / runtime[0:N] - ngrdcol[0:N] / T_cpu[0:N] )**2 )
    rms_error = np.sqrt( np.mean( ( ngrdcol / runtime - ngrdcol / T_cpu )**2 ) )

    return rms_error

def model_throughput(df, column_name):

    ngrdcol = df["ngrdcol"].values  # Access ngrdcol if needed
    runtime = df[column_name].values  # Convert selected column to numpy array


    initial_guess = [2500, 0.002, 1.2]
    result = minimize(objective, initial_guess, args=(ngrdcol, runtime), method='Nelder-Mead')
    #print(result)
    o_tuned, k_tuned, c_tuned = result.x

    print(f"Tuned parameters: o={o_tuned}, k={k_tuned}, c={c_tuned}")
    print(f"Final RMS Error: {result.fun}")

    params = [ o_tuned, k_tuned, c_tuned ]
    #params = [ 1756, .002, 1.13 ]
    T_cpu = model_cpu_time(ngrdcol, runtime, params)


    # L = len(runtime)

    # N_cores = 128
    # N_vsize = 256
    # N_prec  = 64

    # d = 0
    # b = runtime[d] - ( runtime[d+1] - runtime[d] ) / ( ngrdcol[d+1] - ngrdcol[d] ) * ngrdcol[d]

    # T_s = runtime[np.where(ngrdcol == N_cores)] - b
    # T_v = runtime[np.where(ngrdcol == N_vsize*N_cores/N_prec)] - b

    # N_s = np.ceil( ngrdcol / N_cores ) % ( N_vsize / N_prec )
    # N_v = np.floor( ngrdcol / ( N_cores * N_vsize / N_prec ) )

    # f = b + T_s * N_s + T_v * N_v

    # T_cpu = runtime / f


    return pd.Series(T_cpu, index=df.index)


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

def launch_dash_app(grouped_files, shared_variables):
    data = {case: {filename: pd.read_csv(filepath, comment="#") for filename, filepath in files.items()} for case, files in grouped_files.items()}
    x_min = min(df["ngrdcol"].min() for case_data in data.values() for df in case_data.values())
    x_max = max(df["ngrdcol"].max() for case_data in data.values() for df in case_data.values())

    app = Dash(__name__)
    app.title = "Dynamic Plotter"

    plot_style = { "margin-bottom": "2%",
                   "border": "1px solid black",
                   "padding": "1%",
                   "width": "calc(100% - 35px)",
                   "height": "calc(33% - 20px)",
                   "align-items": "center" }

    app.layout = html.Div([
        html.Div([
            html.Div(dcc.Graph(id="plot-columns-per-second"), style=plot_style),
            html.Div(dcc.Graph(id="plot-raw"), style=plot_style),
            html.Div(dcc.Graph(id="fit-plot"), style=plot_style),
            html.Div([
                dcc.Graph(id="plot-custom"),
                dcc.Input(
                    id="custom-function",
                    type="text",
                    placeholder="Enter function (e.g., baseline + 5 * iter_total)",
                    debounce=True,
                    style={"width": "30%", "margin-top": "10px", "text-align": "center"}
                ),
            ], style=plot_style),
        ], style={"flex": "1", "padding-right": "10px", "align-items": "center"}),

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
                            # Default select the first two files in the first group only
                            value=[
                                f"{case}/{filename}"
                                for filename in list(files.keys())[:2]
                            ] if i == 0 and files else []
                        )
                    ],
                    open=(i == 0)  # Default expand the first group only
                )
                for i, (case, files) in enumerate(grouped_files.items())
            ], style={"margin-bottom": "20px"}),
            html.H4("Select variable to plot:"),
            dcc.Dropdown(
                id="variable-dropdown",
                options=[{"label": var, "value": var} for var in shared_variables],
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
            html.Button(
                "Fit",
                id="fit-function-button",
                n_clicks=0,
                style={"margin-top": "10px", "padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
            ),
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
    ], style={"display": "flex"})

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
                temp_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                temp_df["Source"] = f"{case}/{filename}"
                combined_df = pd.concat([combined_df, temp_df])
        fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Source")
        fig.update_layout( xaxis=dict(type=xaxis_scale))
        fig.update_layout( yaxis=dict(type=yaxis_scale))
        return plot_with_enhancements(fig, f"{selected_variable} vs. Number of Grid Columns (Raw Values)")

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
                temp_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                temp_df["Source"] = f"{case}/{filename}"
                temp_df["Columns per Second"] = temp_df["ngrdcol"] / temp_df[selected_variable]
                combined_df = pd.concat([combined_df, temp_df])
        fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Source")
        fig.update_layout( xaxis=dict(type=xaxis_scale))
        fig.update_layout( yaxis=dict(type=yaxis_scale))
        return plot_with_enhancements(fig, "Columns per Second vs. Number of Grid Columns")

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
        fig.update_layout( xaxis=dict(type=xaxis_scale))
        fig.update_layout( yaxis=dict(type=yaxis_scale))
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
    Output("fit-plot", "figure"),
        [
            Input("fit-function-button", "n_clicks"),
            Input({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value")  
        ]
    )
    def update_fit_plot(n_clicks, selected_files, selected_variable, xaxis_scale, yaxis_scale):

        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()

        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:

                temp_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                temp_df["Name"] = f"{case}/{filename}"
                temp_df["Columns per Second"] = temp_df["ngrdcol"] / temp_df[selected_variable]
                combined_df = pd.concat([combined_df, temp_df])

                # Apply model_throughput to the selected variable and add it as "_dup"
                temp_df_dup = temp_df.copy()
                temp_df_dup[selected_variable] = model_throughput(temp_df, selected_variable)
                temp_df_dup["Columns per Second"] = temp_df_dup["ngrdcol"] / temp_df_dup[selected_variable]
                temp_df_dup["Name"] = f"{case}/{filename}_model"
                combined_df = pd.concat([combined_df, temp_df_dup])

        fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Name")
        #fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Name")

        fig.update_layout( xaxis=dict(type=xaxis_scale))
        fig.update_layout( yaxis=dict(type=yaxis_scale))
        return plot_with_enhancements(fig, f"{selected_variable} vs. Number of Grid Columns (Raw Values)")

    app.run_server(debug=True,port=8051)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot shared variables across multiple CSV files.")
    parser.add_argument("-dir", required=True, help="Directory to look for .csv files.")
    args = parser.parse_args()

    csv_files = glob(os.path.join(args.dir, "*.csv"))
    grouped_files = group_files_by_case(csv_files)

    shared_variables, variable_sets = get_shared_variables(
        [file for files in grouped_files.values() for file in files.values()]
    )

    if not shared_variables:
        print("No shared variables found across the CSV files.")
        print("Variable presence per file:")

        # Get all unique variables across files
        all_variables = set.union(*variable_sets.values()) if variable_sets else set()

        # Report missing variables for each file
        for file, variables in variable_sets.items():
            missing_vars = all_variables - variables
            if missing_vars:
                print(f"{file} is missing the variables: {', '.join(sorted(missing_vars))}")

        exit(1)

    launch_dash_app(grouped_files, shared_variables)
