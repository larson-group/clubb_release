import re
import sys
import os
from collections import defaultdict
import numpy as np
import dash
from dash import dcc, html, Input, Output, State
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from joblib import Parallel, delayed
import argparse


nz_max = 134
NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
NZ_RE = re.compile(r"nz_max\s*=\s*(\d+)", re.IGNORECASE)
ITER_RE = re.compile(r"iteration\s*=\s*(\d+)", re.IGNORECASE)


def clean_text(s: str) -> str:
    return re.sub(r"[\x00-\x1F]+", "", s).strip()


def parse_log_file(path: str, debug: bool = True):
    with open(path, "r", errors="replace") as f:
        raw = f.read()
    lines = [line.replace("\x00", "") for line in raw.splitlines()]

    # --- extract nz_max once from whole file ---
    nz_max_val = None
    m = NZ_RE.search(raw)
    if m:
        try:
            nz_max_val = int(m.group(1))
        except ValueError:
            pass

    start_idxs = [i for i, line in enumerate(lines) if "filling holes of:" in line.lower()]
    if not start_idxs:
        return [], nz_max_val
    start_idxs.append(len(lines))

    data = []
    occ_counter = defaultdict(int)
    current_iter = 0  # default before first print

    for bi in range(len(start_idxs) - 1):
        s_idx, e_idx = start_idxs[bi], start_idxs[bi + 1]
        block = lines[s_idx:e_idx]

        # NAME
        name_line = block[0]
        name_part = name_line.split(":", 1)[1] if ":" in name_line else ""
        name = clean_text(name_part) or "(unknown)"

        # THRESHOLD
        threshold = float("nan")
        thresh_line = next((ln for ln in block if "field threshold:" in ln.lower()), None)
        if thresh_line:
            after = thresh_line.lower().split("field threshold:", 1)[1]
            m = NUM_RE.search(after)
            if m:
                try:
                    threshold = float(m.group(0))
                except ValueError:
                    pass

        # FIELD (before)
        field_vals = []
        field_start_idx = None
        for j, ln in enumerate(block):
            if "field:" in ln.lower() and "field after:" not in ln.lower():
                field_start_idx = j
                tail = ln.lower().split("field:", 1)[1]
                field_vals.extend(float(x) for x in NUM_RE.findall(tail))
                break
        if field_start_idx is not None:
            for ln in block[field_start_idx + 1:]:
                if "-- over" in ln.lower():
                    break
                nums = NUM_RE.findall(ln)
                if nums:
                    field_vals.extend(float(x) for x in nums)
        arr = np.array(field_vals, dtype=float)

        # FIELD AFTER
        after_vals = []
        after_start_idx = None
        for j, ln in enumerate(block):
            if "field after:" in ln.lower():
                after_start_idx = j
                tail = ln.lower().split("field after:", 1)[1]
                after_vals.extend(float(x) for x in NUM_RE.findall(tail))
                break
        if after_start_idx is not None:
            for ln in block[after_start_idx + 1:]:
                if "-- over" in ln.lower():
                    break
                nums = NUM_RE.findall(ln)
                if nums:
                    after_vals.extend(float(x) for x in nums)
        arr_after = np.array(after_vals, dtype=float)

        occ = occ_counter[name]
        occ_counter[name] += 1

        data.append({
            "block_index": bi,
            "occurrence": occ,
            "iteration": current_iter,   # << NEW main axis
            "name": name,
            "threshold": threshold,
            "field": arr,
            "field_after": arr_after,
        })

        # --- update iteration if a print is inside this block ---
        for ln in block:
            m = ITER_RE.search(ln)
            if m:
                try:
                    current_iter = int(m.group(1))
                except ValueError:
                    pass
                break

    return data, nz_max_val


# ---------- Build a single figure dict from precomputed record ----------
def build_field_figure_dict(pre):
    fig = go.Figure()

    if pre["field"].size == 0:
        return go.Figure().to_plotly_json()

    # BEFORE (circles, colored by threshold)
    fig.add_trace(go.Scatter(
        x=pre["x"], y=pre["field"],
        mode="markers",
        marker=pre["marker"],
        name="field"
    ))

    # AFTER (triangles at indices that changed)
    if pre["diffs"] and pre["diff_idx"] is not None and len(pre["diff_idx"]):
        fig.add_trace(go.Scatter(
            x=pre["diff_idx"], y=pre["modified_levels"],
            mode="markers",
            marker=pre["after_marker"],
            name="after (diff)"
        ))

    thr = pre["threshold"]
    if np.isfinite(thr):
        fig.add_hline(y=thr, line_dash="dot", line_color="black")

    fig.update_yaxes(type="log")
    fig.update_layout(
        title=f"{pre['name']} | iter {pre['iteration']} | thr={thr}",
        height=500
    )
    return fig.to_plotly_json()


# ---------------------- MAIN ----------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run fill_holes app with port and input file."
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8050,
        help="Port number for the Dash app (default: 8050)"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the text file (e.g., arm_fill)"
    )

    args = parser.parse_args()

    filepath = args.input_file
    data, nz_max_val = parse_log_file(filepath, debug=False)
    if not data:
        print("[ERROR] No data parsed from file.")
        sys.exit(1)

    dash_port = args.port
    nz_max = nz_max_val if nz_max_val is not None else 134

    # Build indices and precompute lightweight arrays/masks once
    by_name_iter = {}
    names_set = set()

    for rec in data:
        arr, arr_after, thr = rec["field"], rec["field_after"], rec["threshold"]
        name, iteration = rec["name"], rec["iteration"]

        x = np.arange(len(arr))
        labels = np.where(arr < thr, "below", "above")
        marker = dict(
            symbol="circle",
            size=10,
            color=["red" if l == "below" else "blue" for l in labels]
        )

        diffs = False
        diff_idx = None
        modified_levels = None
        after_marker = None
        if len(arr_after):
            min_len = min(len(arr), len(arr_after))
            diff_idx = np.where(arr_after[:min_len] != arr[:min_len])[0]
            if len(diff_idx):
                diffs = True
                modified_levels = arr_after[diff_idx]
                after_colors = np.where(modified_levels < thr, "red", "green")
                after_marker = dict(symbol="triangle-up", size=12, color=after_colors)

        by_name_iter[(name, iteration)] = {
            "name": name,
            "iteration": iteration,
            "field": arr,
            "field_after": arr_after,
            "threshold": thr,
            "x": x,
            "marker": marker,
            "diffs": diffs,
            "diff_idx": diff_idx,
            "modified_levels": modified_levels,
            "after_marker": after_marker,
        }

        names_set.add(name)

    # Summary data
    summary_records = []
    for rec in data:
        below = int(np.sum(rec["field"] < rec["threshold"]))
        summary_records.append({
            "name": rec["name"],
            "iteration": rec["iteration"],
            "below_count": below,
        })
    summary_df = pd.DataFrame(summary_records)

    # Build first-plot figures in parallel
    by_name_keys = defaultdict(list)
    for (nm, it) in by_name_iter.keys():
        by_name_keys[nm].append((nm, it))
    for nm in by_name_keys:
        by_name_keys[nm].sort(key=lambda t: t[1])

    def build_all_for_field(name, keys):
        out = []
        for key in keys:
            fig_dict = build_field_figure_dict(by_name_iter[key])
            out.append((key, fig_dict))
        return out

    num_jobs = min(len(by_name_keys), max(1, os.cpu_count() or 1))
    results = Parallel(n_jobs=num_jobs)(
        delayed(build_all_for_field)(name, keys)
        for name, keys in by_name_keys.items()
    )

    prebuilt_figs = {}
    for chunk in results:
        for key, fig_dict in chunk:
            prebuilt_figs[key] = fig_dict
    by_name_iter.clear()

    # ---------- Dash app ----------
    names = sorted(names_set)
    default_name = names[0]

    global_max_iter = max(rec["iteration"] for rec in data)
    default_max = global_max_iter
    default_iter = 0

    app = dash.Dash(__name__)
    app.layout = html.Div(
        [
            html.H1("Simulation Field Visualizer"),

            html.Div([
                html.Label("Field NAME"),
                dcc.Dropdown(
                    id="field-name",
                    options=[{"label": n, "value": n} for n in names],
                    value=default_name,
                    clearable=False,
                ),
            ], style={"maxWidth": "480px"}),

            html.Br(),

            html.Div([
                html.Label("Iteration"),
                dcc.Slider(
                    id="iteration-slider",
                    min=0,
                    max=default_max,
                    step=1,
                    value=default_iter,
                    marks=None,
                    tooltip={"placement": "bottom", "always_visible": False},
                    updatemode="drag",
                ),
            ], style={"width": "100%"}),

            html.Div([
                html.Button("Play/Pause", id="play-button", n_clicks=0),
                html.Label("Playback speed (ms per step):"),
                dcc.Slider(
                    min=50, max=2000, step=100, value=500,
                    id="speed-slider",
                    tooltip={"placement": "bottom"}
                ),
                dcc.Interval(id="play-interval", interval=500, n_intervals=0, disabled=True)
            ], style={"width": "30%", "padding": "20px"}),

            dcc.Graph(id="field-plot"),
            dcc.Graph(id="below-mask-plot"),
            dcc.Graph(id="summary-plot"),
        ]
    )

    # Playback toggle
    @app.callback(
        Output("play-interval", "disabled"),
        Input("play-button", "n_clicks"),
        State("play-interval", "disabled")
    )
    def toggle_play(n_clicks, disabled):
        if n_clicks is None:
            return True
        return not disabled

    # Update playback speed
    @app.callback(
        Output("play-interval", "interval"),
        Input("speed-slider", "value"),
    )
    def update_speed(speed_val):
        return speed_val

    # Increment slider on interval
    @app.callback(
        Output("iteration-slider", "value", allow_duplicate=True),
        Input("play-interval", "n_intervals"),
        State("iteration-slider", "value"),
        State("iteration-slider", "max"),
        prevent_initial_call=True,
    )
    def advance_iteration(n, current, max_val):
        if current is None:
            return 0
        return (current + 1) % (max_val + 1)

    # First plot
    @app.callback(
        Output("field-plot", "figure"),
        Input("field-name", "value"),
        Input("iteration-slider", "value"),
    )
    def update_field_plot(name, iteration):
        fig_dict = prebuilt_figs.get((name, iteration))
        return fig_dict or go.Figure().to_plotly_json()

    # Second plot
    @app.callback(
        Output("below-mask-plot", "figure"),
        Input("field-name", "value"),
    )
    def update_below_mask_plot(name):
        rows = []
        for rec in data:
            if rec["name"] != name:
                continue
            thr = rec["threshold"]
            for idx in np.where(rec["field"] < thr)[0]:
                rows.append({"iteration": rec["iteration"], "field_index": idx, "kind": "initial"})
            for idx in np.where(rec["field_after"] < thr)[0]:
                rows.append({"iteration": rec["iteration"], "field_index": idx, "kind": "unfixed"})
        if not rows:
            return go.Figure()
        df = pd.DataFrame(rows)
        fig = px.scatter(
            df, x="iteration", y="field_index",
            color="kind",
            color_discrete_map={"initial": "green", "unfixed": "red"},
            title=f"Indices Below Threshold over Time for {name}",
        )
        fig.update_traces(marker=dict(size=7))
        max_iter = max(rec["iteration"] for rec in data if rec["name"] == name)
        fig.update_layout(
            height=500,
            xaxis=dict(range=[-5, max_iter+5]),
            yaxis=dict(range=[-5, nz_max+5])
        )
        return fig

    # Third plot
    @app.callback(
        Output("summary-plot", "figure"),
        Input("field-name", "value")
    )
    def update_summary_plot(_):
        fig = px.line(
            summary_df,
            x="iteration", y="below_count",
            color="name", markers=True,
            title="Count of Values Below Threshold per Iteration (all fields)",
        )
        fig.update_layout(height=500)
        return fig

    app.run(debug=False, port=dash_port)
