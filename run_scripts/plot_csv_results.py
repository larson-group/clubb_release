#!/bin/python3

import argparse
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
import warnings
from itertools import chain

# ======================================== Reused html styles ========================================
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


# ======================================== Cache penalty functions ========================================
def loglog( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return c * np.log(np.log( np.exp(k*(x - o)) + 1 ) + 1 ) + 1

def sigmoid( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return 1 + c / (1 + np.exp(-k * (x - o)))

def sqrtlog( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return c * np.sqrt( np.log( np.exp(k*(x - o)) + 1 ) ) + 1

def logsig( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return  c * np.log(np.log( np.exp(x - o) + 1 ) + 1 ) / (1 + np.exp(-k * ( x - o )) ) + 1

def logcosh( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return c * np.log(np.log(np.cosh(k * np.maximum(x - o, 0.0) )) + 1) + 1

def logcoshsig( c, k, o, x ):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    return c * np.log( np.log(np.cosh(np.maximum(x-o, 0.0) )) / (1 + np.exp(-k * ( x - o )) ) + 1) + 1

#cache_pen_funcs = [ loglog, logcosh, sigmoid, sqrtlog ]
cache_pen_funcs = {
    "loglog": loglog,
    "logcosh": logcosh,
    "logsig": logsig,
    "sigmoid": sigmoid,
    "logcoshsig": logcoshsig,
    "sqrtlog": sqrtlog
}

def rms_error( ngrdcol, runtime, T_model, mode="abs_percent_cps" ):

    abs_diff_time = np.abs( runtime - T_model )
    abs_percent_time = 100 * ( runtime - T_model ) / runtime

    cps = ngrdcol / runtime
    cps_model = ngrdcol / T_model

    abs_percent_diff_cps = 100 * ( cps - cps_model ) / cps
    abs_diff_cps = np.abs( cps - cps_model )

    if mode == "abs_percent_cps":
        rms_error = np.sqrt( np.mean( abs_percent_diff_cps**2 ) )
    elif mode == "abs_cps":
        rms_error = np.sqrt( np.mean( abs_diff_cps**2 ) )


    return rms_error



# ======================================== VCPU ========================================

vcpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "T_v", "id": "T_v_val"},
    {"name": "T_r", "id": "T_r_val"},
    {"name": "b", "id": "b_val"},
    #{"name": "b_est", "id": "b_est_val"},
    {"name": "c", "id": "c_val"},
    {"name": "k", "id": "k_val"},
    {"name": "o", "id": "o_val"},
    {"name": "Cache Pen Func", "id": "cp_func"},
    {"name": "RMSE", "id": "rms_error"},
]


def vcpu_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):
    _, rms_error = model_vcpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )
    return rms_error
def model_vcpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    # Model time
    cols_per_core = ngrdcol / N_tasks
    flops_per_vop = N_vsize / N_prec

    T_r, T_v, b, c, k, o = param_scale * params

    #T_r = runtime[np.where(cols_per_core == 1)] - b
    #T_v = runtime[np.where(cols_per_core == flops_per_vop)] - b

    N_s = np.ceil( cols_per_core ) % ( flops_per_vop )
    N_v = np.floor( cols_per_core / ( flops_per_vop ) )

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20

    T_vcpu = ( T_r * N_s + T_v * N_v ) * cp_func( c, k, o, avg_array_size_MB ) + b

    return T_vcpu, rms_error( ngrdcol, runtime, T_vcpu )


# ======================================== CPU ========================================

cpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m", "id": "m_val"},
    #{"name": "m_est", "id": "m_est_val"},
    {"name": "b", "id": "b_val"},
    #{"name": "b_est", "id": "b_est_val"},
    {"name": "c", "id": "c_val"},
    {"name": "k", "id": "k_val"},
    {"name": "o", "id": "o_val"},
    {"name": "Cache Pen Func", "id": "cp_func"},
    {"name": "RMSE", "id": "rms_error"},
]

def cpu_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):
    _, rms_error = model_cpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )
    return rms_error
def model_cpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    # Model time
    m, b, c, k, o = param_scale * params

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20
 
    T_cpu =  ( m * ngrdcol + b ) * cp_func( c, k, o, avg_array_size_MB )

    return T_cpu, rms_error( ngrdcol, runtime, T_cpu )

# ======================================== CPU Batched ========================================

cpu_batched_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m_ik", "id": "m_ik_val"},
    {"name": "m_k", "id": "m_k_val"},
    #{"name": "m_est", "id": "m_est_val"},
    {"name": "b", "id": "b_val"},
    #{"name": "b_est", "id": "b_est_val"},
    {"name": "c", "id": "c_val"},
    {"name": "k", "id": "k_val"},
    {"name": "o", "id": "o_val"},
    {"name": "Cache Pen Func", "id": "cp_func"},
    {"name": "RMSE", "id": "rms_error"},
]

def cpu_batched_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):
    params = params * param_scale
    _, rms_error = model_cpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )
    return rms_error
def model_cpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    # Model time
    m_ik, m_k, b, c, k, o = params

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20
 
    T_cpu =  ( m_ik * ngrdcol * N_vlevs + m_k * N_vlevs + b ) * cp_func( c, k, o, avg_array_size_MB ) 

    return T_cpu, rms_error( ngrdcol, runtime, T_cpu )


# ======================================== GPU ========================================

gpu_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m", "id": "m_val"},
    {"name": "m_est", "id": "m_est_val"},
    {"name": "b", "id": "b_val"},
    {"name": "b_est", "id": "b_est_val"},
    {"name": "RMSE", "id": "rms_error"},
]

def gpu_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs ):
    _, rms_error = model_gpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs )
    return rms_error
def model_gpu_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs):

    m, b = param_scale * params

    T_gpu = b + m * ngrdcol

    return T_gpu, rms_error( ngrdcol, runtime, T_gpu, mode="abs_cps" )


# ======================================== Batched GPU ========================================

gpu_batched_model_columns = [
    {"name": "Name", "id": "name"},
    {"name": "m_ik", "id": "m_ik_val"},
    {"name": "m_k", "id": "m_k_val"},
    {"name": "b", "id": "b_val"},
    #{"name": "b_est", "id": "b_est_val"},
    {"name": "d", "id": "d_val"},
    #{"name": "m_i", "id": "m_i_val"},
    {"name": "RMSE", "id": "rms_error"},
]

def gpu_batched_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs ):
    _, rms_error = model_gpu_batched_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs )
    return rms_error
def model_gpu_batched_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs):

    m_ik, m_k, b = param_scale * params

    T_bgpu = m_ik * ngrdcol * N_vlevs + m_k * N_vlevs + b

    return T_bgpu, rms_error( ngrdcol, runtime, T_bgpu, mode="abs_cps" )


def model_throughputs(ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, model_version, selected_cp_funcs=None):


    #============================== VCPU Model ==============================

    if model_version == "vcpu":

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]
        T_r_est = runtime[np.where(cols_per_core == 1)] - b_est
        T_v_est = runtime[np.where(cols_per_core == flops_per_vop)] - b_est

        initial_guess = [ T_r_est[0], T_v_est[0], b_est, 1.0, 1.0, 1.0 ]

        print(initial_guess)

        param_scale = np.abs( initial_guess )
        bounds = [ (None,None) ] * len(initial_guess)

        cache_pen_best = None
        rms_min = None
        params_opt = None

        for name in selected_cp_funcs:

            cp_func = cache_pen_funcs[name]

            print(f" - trying cache pen: {cp_func.__name__}")

            result = minimize(  vcpu_objective, initial_guess, 
                                args=(param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func), 
                                #method='BFGS'
                                method='L-BFGS-B',
                                bounds = bounds)

            rms_error = result.fun
            print(f" -- rms_error = {rms_error}")

            if rms_min is None or rms_error < rms_min:
                cache_pen_best = cp_func
                rms_min = rms_error
                params_opt = result.x


        T_cpu, _ = model_vcpu_time(params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        params_opt = param_scale * params_opt

        fit_params = {
            "T_v_val": params_opt[0], 
            "T_r_val": params_opt[1], 
            "b_val": params_opt[2],  
            #"b_est_val": b_est, 
            "c_val": params_opt[3],
            "k_val": params_opt[4],
            "o_val": params_opt[5], 
            "cp_func": cache_pen_best.__name__, 
            "rms_error" : rms_error
        }

    #============================== CPU Model ==============================

    if model_version == "cpu":

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec
        first_vpoint = np.where(cols_per_core == flops_per_vop)[0][0]

        m_est = ( runtime[first_vpoint] - b_est ) / ( ngrdcol[first_vpoint] )

        initial_guess = [ m_est, b_est, 1.0, 1.0, 1.0 ]

        param_scale = np.abs( initial_guess )
        bounds = [ (None,None) ] * len(initial_guess)

        cache_pen_best = None
        rms_min = None
        params_opt = None

        for name in selected_cp_funcs:

            cp_func = cache_pen_funcs[name]

            print(f" - trying cache pen: {cp_func.__name__}")


            result = minimize(
                cpu_objective,
                initial_guess,
                args=(param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func),
                #method='BFGS'
                method='L-BFGS-B',
                bounds = bounds
            )

            rms_error = result.fun
            print(f" -- rms_error = {rms_error}")

            if rms_min is None or rms_error < rms_min:
                cache_pen_best = cp_func
                rms_min = rms_error
                params_opt = result.x


        T_cpu, _ = model_cpu_time(params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

        params_opt = param_scale * params_opt

        fit_params = {
            "m_val": params_opt[0],
            "m_est_val": m_est, 
            "b_val": params_opt[1], 
            "b_est_val": b_est, 
            "c_val": params_opt[2], 
            "k_val": params_opt[3], 
            "o_val": params_opt[4], 
            "cp_func": cache_pen_best.__name__, 
            "rms_error" : rms_min
        }

    #============================== GPU Model ==============================

    if model_version == "gpu":

        b_est = runtime[-2] - ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] ) * ngrdcol[-2]
        m_est = ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] )

        initial_guess = [ m_est, b_est ]

        param_scale = np.abs( initial_guess )
        bounds = [ (0,None) ] * len(initial_guess)

        result = minimize(  gpu_objective, 
                            initial_guess, 
                            args=(param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs), 
                            #method='Nelder-Mead' )
                            method='L-BFGS-B',
                            bounds = bounds )

        rms_error = result.fun
        params_opt = result.x

        T_cpu, _ = model_gpu_time(param_scale, params_opt, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs)

        params_opt = param_scale * params_opt

        fit_params = {
            "m_val": params_opt[0], 
            "m_est_val": m_est, 
            "b_val": params_opt[1], 
            "b_est_val": b_est, 
            "rms_error" : rms_error
        }

    #============================== GPU Batched Model ==============================

    if model_version == "gpu_batched":

        b_est = runtime[-2] - ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] ) * ngrdcol[-2]
        m_ik_est = ( runtime[-1] - runtime[-2] ) / ( ngrdcol[-1] - ngrdcol[-2] ) * ( 1 / N_vlevs[-1])
        m_k_est = b_est / N_vlevs[-1]
        #d_est = 100000.0
        #m_i_est = m_ik_est

        initial_guess = [ m_ik_est, m_k_est, b_est ]

        param_scale = np.abs( initial_guess )
        bounds = [ (0,None) ] * len(initial_guess)

        params_opt = None

        result = minimize(  gpu_batched_objective, 
                            initial_guess, 
                            args=( param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs ), 
                            #method='Nelder-Mead' )
                            method='L-BFGS-B',
                            bounds = bounds )

        rms_error = result.fun
        params_opt = result.x

        T_cpu, _ = model_gpu_batched_time( params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs )

        params_opt = param_scale * params_opt

        fit_params = {
            "m_ik_val": params_opt[0], 
            "m_k_val": params_opt[1], 
            "b_val": params_opt[2], 
            "b_est_val": b_est, 
            #"d_val": params_opt[3], 
            #"m_i_val": params_opt[3], 
            "rms_error" : rms_error
        }

    #============================== CPU Batched Model ==============================

    if model_version == "cpu_batched":

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec
        first_vpoint = np.where(cols_per_core == flops_per_vop)[0][0]

        m_ik_est = ( runtime[first_vpoint] - b_est ) / ( ngrdcol[first_vpoint] ) * ( 1 / N_vlevs[-1])
        m_k_est = b_est / N_vlevs[-1]

        param_scale = np.abs( [ m_ik_est, m_k_est, b_est, 1.0, 1.0, 1.0 ] )

        initial_guess = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]

        bounds = [  (0,None),       # m_ik
                    (0,None),       # m_k
                    (0,None),       # b
                    (0,None),       # c
                    (0,None),       # k
                    (None,None)     # o
                ]

        cache_pen_best = None
        rms_min = None
        params_opt = None

        for name in selected_cp_funcs:

            cp_func = cache_pen_funcs[name]

            print(f" - trying cache pen: {cp_func.__name__}")

            result = minimize(  cpu_batched_objective, 
                                initial_guess, 
                                args=( param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ), 
                                #method='Nelder-Mead' )
                                method='L-BFGS-B',
                                bounds = bounds )

            rms_error = result.fun
            scaled_params = param_scale * result.x
            formatted_params = " ".join([f"{p:.3e}" for p in scaled_params])
            print(f" -- rms_error = {rms_error} -- {formatted_params}")
            #print(f" -- rms_error = {rms_error} -- {(param_scale *result.x):.3e}")

            if rms_min is None or rms_error < rms_min:
                cache_pen_best = cp_func
                rms_min = rms_error
                params_opt = result.x

        params_opt = param_scale * params_opt

        T_cpu, _ = model_cpu_batched_time( params_opt, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best )



        # Bootstrap procedure
        params_bootstrap = []

        n_points = len(runtime)
        rng = np.random.default_rng()

        for _ in range(100):
            indices = rng.integers(0, n_points, size=n_points)  # resample with replacement
            ngrdcol_bs = ngrdcol[indices]
            runtime_bs = runtime[indices]
            N_prec_bs = N_prec[indices]
            N_vlevs_bs = N_vlevs[indices]
            #N_tasks_bs = N_tasks if np.ndim(N_tasks) == 0 else N_tasks[indices]
            #N_vsize_bs = N_vsize if np.ndim(N_vsize) == 0 else N_vsize[indices]

            try:
                result_bs = minimize(
                    cpu_batched_objective, 
                    initial_guess, 
                    args=(param_scale, ngrdcol_bs, runtime_bs, N_tasks, N_vsize, N_prec_bs, N_vlevs_bs, cache_pen_best),
                    method='L-BFGS-B',
                    bounds=bounds
                )
                params_bootstrap.append(result_bs.x * param_scale)  # save rescaled params
            except Exception as e:
                # Skip failed fits
                print(f"Bootstrap fit failed: {e}")
                continue

        params_bootstrap = np.array(params_bootstrap)

        confidence_level = 90

        # Now get 5th and 95th percentile (for 90% confidence interval)
        params_low = np.percentile(params_bootstrap, (100 - confidence_level) / 2, axis=0)
        params_high = np.percentile(params_bootstrap, 100 - (100 - confidence_level) / 2, axis=0)

        # Assemble fit_params dictionary
        fit_params = {
            "m_ik_val": f"{params_opt[0]:.2e} \n({params_low[0]:.2e}, {params_high[0]:.2e})",
            "m_k_val": f"{params_opt[1]:.2e} \n({params_low[1]:.2e}, {params_high[1]:.2e})",
            "b_val": f"{params_opt[2]:.2e} \n({params_low[2]:.2e}, {params_high[2]:.2e})",
            "c_val": f"{params_opt[3]:.2e} \n({params_low[3]:.2e}, {params_high[3]:.2e})",
            "k_val": f"{params_opt[4]:.2e} \n({params_low[4]:.2e}, {params_high[4]:.2e})",
            "o_val": f"{params_opt[5]:.2e} \n({params_low[5]:.2e}, {params_high[5]:.2e})",
            "cp_func": cache_pen_best.__name__,
            "rms_error": rms_min
        }

    return T_cpu, fit_params


def launch_dash_app(dir_name, grouped_files, all_variables):

    data = {case: {filename: pd.read_csv(filepath, comment="#") for filename, filepath in files.items()} for case, files in grouped_files.items()}

    #app = Dash(__name__)
    app = Dash(
            __name__,
            requests_pathname_prefix='/plots/',
            routes_pathname_prefix='/plots/'
    )

    app.title = "Dynamic Plotter"
                
    app.layout = html.Div([

        # ======================================== Plots ========================================
        html.Div([

            # Columns per second plot
            html.Div(dcc.Graph(id="plot-columns-per-second", style=graph_style), style=plot_div_style),

            # Time plot
            html.Div(dcc.Graph(id="plot-raw", style=graph_style), style=plot_div_style),

            # Fit plot
            html.Div([

                # Plot
                dcc.Graph(id="fit-plot", style=graph_style), 

                # VCPU table
                html.Div([
                    html.Label("Vector CPU Model"),
                    dash_table.DataTable(
                        id='vcpu-param-table',
                        columns=vcpu_model_columns,
                        style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                        data=[],  # Initial empty table
                        editable=True,
                        row_deletable=False  # Allows deleting rows directly
                    ),
                ], style={"margin-bottom": "20px"}),

                # CPU table
                html.Div([
                    html.Label("Simple CPU Model"),
                    dash_table.DataTable(
                        id='cpu-param-table',
                        columns=cpu_model_columns,
                        style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                        data=[],  # Initial empty table
                        editable=True,
                        row_deletable=False  # Allows deleting rows directly
                    ),
                ], style={"margin-bottom": "20px"}),
                
                # GPU table
                html.Div([
                    html.Label("GPU Model"),
                    dash_table.DataTable(
                        id='gpu-param-table',
                        columns=gpu_model_columns,
                        style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                        data=[],  # Initial empty table
                        editable=True,
                        row_deletable=False  # Allows deleting rows directly
                    ),
                ], style={"margin-bottom": "20px"}),

            ], style=fit_plot_div_style),

            # Batched fit plot
            html.Div([

                # Plot
                dcc.Graph(id="fit-plot-batched", style=graph_style), 

                # GPU Batched Table
                html.Label("GPU Model"),
                dash_table.DataTable(
                    id='gpu-batched-param-table',
                    columns=gpu_batched_model_columns,
                    style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                    data=[],  # Initial empty table
                    editable=True,
                    row_deletable=False  # Allows deleting rows directly
                ),

                # CPU Batched Table
                html.Label("CPU Model"),
                dash_table.DataTable(
                    id='cpu-batched-param-table',
                    columns=cpu_batched_model_columns,
                    style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                    data=[],  # Initial empty table
                    editable=True,
                    row_deletable=False  # Allows deleting rows directly
                ),

            ], style=fit_plot_div_style),

            # Custom function plot
            html.Div([

                # Plot
                dcc.Graph(id="plot-custom", style=graph_style),

                # Custom function input
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
                    "height": "fit-content" }
        ),


        # ======================================== Controls ========================================
        html.Div([

            # File selection
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

            # Variable selection
            html.H4("Select variable to plot:"),
            dcc.Dropdown(
                id="variable-dropdown",
                options=[{"label": var, "value": var} for var in all_variables],
                #value = "mainloop_r" if "_gptl" in dir_name else "compute_i",],
                value = "compute_i",
                style={"margin-bottom": "20px"}
            ),

            # Clipboard copy 
            html.Button(
                "Copy CSV Data to Clipboard",
                id="copy-csv-button",
                n_clicks=0,
                style={"margin-top": "10px", "padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
            ),

            # x-axis scale 
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

            # y-axis scale 
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

            # Fit function controls
            html.Div([
                
                html.Div([
                    html.Label("Fit Function Variables:", style={"font-weight": "bold"}),
                ], style={"margin-bottom": "8px"}),

                html.Div([

                    # Time / cps choice
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
                    ], style={"margin-bottom": "8px"}),

                    # Cache pen func choice
                    html.Div([
                        html.Details(
                            [
                                html.Summary("cache penalty functions"),
                                dcc.Checklist(
                                    id={"type": "cp_func_checkbox"},
                                    options=[
                                        {"label": name, "value": name}
                                        for name in cache_pen_funcs
                                    ],
                                    value=["logcosh"],
                                    labelStyle={"display": "inline-block", "width": "25%"}  # 2 columns
                                )
                            ],
                            open=True
                        )
                    ], style={"margin-bottom": "8px"}),


                    html.Div([
                        html.Label("Vector Length: ", style={"margin-top": "5px"}),
                        dcc.Input(id="N_vsize", type="number", value=256, style={"width": "20%"}),
                    ], style={"margin-bottom": "8px"}),

                    # Fit button
                    html.Button(
                        "Fit",
                        id="fit-function-button",
                        n_clicks=0,
                        style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
                    ),

                    # Batched fit button
                    html.Button(
                        "Fit Batched",
                        id="fit-batched-function-button",
                        n_clicks=0,
                        style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
                    ),

                ], style={"display": "flex", "flex-direction": "column", "gap": "5px", "margin-bottom": "10px"}),


            ], style={"margin-top": "10px", "padding": "10px", "border": "1px solid #000", "background-color": "#f9f9f9"}),

        ], style={
            "width": "20%",
            "padding": "10px",
            "border": "1px solid black",
            "position": "sticky",
            "top": "0",
            "height": "fit-content",
            "background-color": "white",
        })
    ], style={"display": "flex", "flexDirection": "row", "width": "100%", "height": "auto"})


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

    # ======================================== Time plot ========================================
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

        # Flatten the list of lists into a single list of selected filenames
        selected_flat = list(chain.from_iterable(selected_files))

        if not selected_flat:
            raise PreventUpdate

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


    # ======================================== CPS plot ========================================
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

        # Flatten the list of lists into a single list of selected filenames
        selected_flat = list(chain.from_iterable(selected_files))

        if not selected_flat:
            raise PreventUpdate

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


    # ======================================== Custom function plot ========================================

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
    
    def evaluate_function(func, df):
        try:
            local_vars = {col: df[col] for col in df.columns}
            return eval(func, {"__builtins__": None}, {**local_vars, **safe_functions})
        except Exception as e:
            print(f"Error evaluating function '{func}': {e}")
            return None
            
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

        # Flatten the list of lists into a single list of selected filenames
        selected_flat = list(chain.from_iterable(selected_files))

        if not selected_flat:
            raise PreventUpdate

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


    # ======================================== Copy button function ========================================
    @app.callback(
    Output("copy-csv-button", "children"),  # Update button text to confirm copy
        [
            Input("copy-csv-button", "n_clicks"),
            State({"type": "file-checkbox", "case": ALL}, "value"),
            State("variable-dropdown", "value")
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



    # ======================================== Fit plot ========================================
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
            State("N_vsize", "value"),
            State({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("fit-plot-mode", "value"),
            State({"type": "cp_func_checkbox"}, "value"),
        ],
    )
    def update_fit_plot(n_clicks, N_vsize, selected_files, selected_variable, xaxis_scale, yaxis_scale, plot_mode, selected_cp_funcs):
        
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
                original_df["Name"] = f"{filename}"
                original_df["Columns per Second"] = original_df["ngrdcol"] / original_df[selected_variable]

                combined_df = pd.concat([combined_df, original_df])
                
                # Get the number of cores used 
                N_tasks = data[case][filename]["tasks"][0]
                N_vlevs = data[case][filename]["nz"][0]

                # Apply model_throughput to the selected variable and add it as "_dup"
                vcpu_df = original_df.copy()
                cpu_df = original_df.copy()
                gpu_model_df = original_df.copy()

                if "_sp_" in filename:
                    N_prec = 32
                else:
                    N_prec = 64
                
                if any(gpu_name in filename for gpu_name in ["A100", "V100", "H100"]):

                    print(f"modeling {filename} with gpu")
                    T_gpu, gpu_params = model_throughputs(  gpu_model_df["ngrdcol"].values, 
                                                            gpu_model_df[selected_variable].values, 
                                                            N_tasks, 
                                                            N_vsize, 
                                                            N_prec, 
                                                            N_vlevs, 
                                                            "gpu")

                    gpu_model_df[selected_variable] = T_gpu
                    gpu_model_df["Columns per Second"] = gpu_model_df["ngrdcol"] / gpu_model_df[selected_variable]
                    gpu_model_df["Name"] = f"{filename}_gmodel"

                    gpu_table_data.append({
                        col["id"]: f"{gpu_params.get(col['id']):.3e}" if isinstance(gpu_params.get(col["id"]), (int, float)) else gpu_params.get(col["id"])
                        for col in gpu_model_columns
                    })
                    gpu_table_data[-1]["name"] = f"{case}/{filename}"

                    combined_df = pd.concat([combined_df, gpu_model_df])

                else:

                    print(f"modeling {filename} with vcpu")
                    T_vcpu, vcpu_params = model_throughputs( vcpu_df["ngrdcol"].values, 
                                                             vcpu_df[selected_variable].values, 
                                                             N_tasks, 
                                                             N_vsize, 
                                                             N_prec, 
                                                             N_vlevs, 
                                                             "vcpu",
                                                             selected_cp_funcs )

                    vcpu_df[selected_variable] = T_vcpu
                    vcpu_df["Columns per Second"] = vcpu_df["ngrdcol"] / vcpu_df[selected_variable]
                    vcpu_df["Name"] = f"{case}/{filename}_vmodel"

                    vcpu_table_data.append({
                        col["id"]: f"{vcpu_params.get(col['id']):.3e}" if isinstance(vcpu_params.get(col["id"]), (int, float)) else vcpu_params.get(col["id"])
                        for col in vcpu_model_columns
                    })
                    vcpu_table_data[-1]["name"] = f"{case}/{filename}"


                    print(f"modeling {filename} with cpu")
                    T_cpu, cpu_params = model_throughputs( cpu_df["ngrdcol"].values, 
                                                           cpu_df[selected_variable].values, 
                                                           N_tasks, 
                                                           N_vsize, 
                                                           N_prec, 
                                                           N_vlevs, 
                                                           "cpu",
                                                           selected_cp_funcs)
                    
                    cpu_df[selected_variable] = T_cpu
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
            fig.update_layout( yaxis_title = "Throughput (columns / second)" )
        else:
            fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Name")
            fig.update_layout( yaxis_title = f"Runtime of {selected_variable} (seconds)" )

        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            xaxis_title = "Batch size (Columns)", 
            autosize=True,
            uirevision='constant'
        )
        fig = plot_with_enhancements(fig, f"{selected_variable} vs. Number of Grid Columns (Raw Values)")

        return fig.to_dict(), vcpu_table_data, cpu_table_data, gpu_table_data


    # ======================================== Batched fit plot ========================================
    @app.callback(
        [
            Output("fit-plot-batched", "figure"),
            #Output('vcpu-param-table', 'data'),
            #Output('cpu-param-table', 'data'),
            Output('gpu-batched-param-table', 'data'),
            Output('cpu-batched-param-table', 'data'),
        ],
        [
            Input("fit-batched-function-button", "n_clicks"),
            #Input("N_tasks", "value"),
            State("N_vsize", "value"),
            State({"type": "file-checkbox", "case": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("fit-plot-mode", "value"),
            State({"type": "cp_func_checkbox"}, "value")
        ],
        prevent_initial_call=True
    )
    def update_fit_plot_batched(n_clicks, N_vsize, selected_files, selected_variable, xaxis_scale, yaxis_scale, plot_mode, selected_cp_funcs):
        
        if n_clicks == 0 or not selected_files:
            raise PreventUpdate

        flat_files = [item.split("/") for sublist in selected_files for item in sublist]
        combined_df = pd.DataFrame()

        gpu_batched_table_data = []
        cpu_batched_table_data = []

        N_tasks = np.empty(0)
        N_vlevs = np.empty(0)
        runtimes = np.empty(0)
        ngrdcols = np.empty(0)
        N_prec = np.empty(0)
        filenames = ""

        gpu_fit = False

        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:
                

                # Get the number of cores used 
                N_tasks  = np.concatenate( [N_tasks, np.array( data[case][filename]["tasks"] )] )
                N_vlevs  = np.concatenate( [N_vlevs, np.array( data[case][filename]["nz"] )] )
                runtimes = np.concatenate( [runtimes, np.array( data[case][filename][selected_variable] )] )
                ngrdcols = np.concatenate( [ngrdcols, np.array( data[case][filename]["ngrdcol"] )] )

                if "_sp_" in filename:
                    N_prec = np.concatenate( [N_prec, np.array([32]*len(data[case][filename]["nz"])) ])
                else:
                    N_prec = np.concatenate( [N_prec, np.array([64]*len(data[case][filename]["nz"])) ])

                # If any of the files in the batch are GPU, just do a GPU fit
                if any(gpu_name in filename for gpu_name in ["A100", "V100", "H100"]):
                    gpu_fit = True 
                
        # This will need changing if we want to batch fit single precision results
        #N_prec = 64

        if gpu_fit:
            T_model, params_opt = model_throughputs(  ngrdcols, 
                                                    runtimes, 
                                                    N_tasks, 
                                                    N_vsize, 
                                                    N_prec, 
                                                    N_vlevs, 
                                                    "gpu_batched")
        else:
            T_model, params_opt = model_throughputs(  ngrdcols, 
                                                    runtimes, 
                                                    N_tasks, 
                                                    N_vsize, 
                                                    N_prec, 
                                                    N_vlevs, 
                                                    "cpu_batched",
                                                    selected_cp_funcs )

        start_idx = 0
        for case, filename in flat_files:
            if case in data and filename in data[case] and selected_variable in data[case][filename].columns:

                original_df = data[case][filename][["ngrdcol", selected_variable]].copy()
                #original_df["Name"] = f"{filename}"
                original_df["Name"] = re.sub(r'(_gptl|_derecho)(?=(_|$))', '', f"{filename}")



                grid_boxes = original_df["ngrdcol"] * data[case][filename]["nz"]
                
                original_df["Grid Boxes"] = grid_boxes
                original_df["Columns per Second"] = original_df["ngrdcol"] / original_df[selected_variable]
                original_df["Grid Boxes per Second"] = grid_boxes / original_df[selected_variable]

                model_df = original_df.copy()

                end_idx = start_idx + len(original_df["ngrdcol"])
                model_df[selected_variable] = T_model[start_idx:end_idx]

                model_df["Grid Boxes per Second"] = grid_boxes / model_df[selected_variable]
                model_df["Columns per Second"] = model_df["ngrdcol"] / model_df[selected_variable]
                model_df["Name"] = original_df["Name"]
                
                original_df["Source"] = "original"
                model_df["Source"] = "model"
                combined_df = pd.concat([combined_df, original_df, model_df])

                # Calculate RMS for specific model
                cps = original_df["ngrdcol"] / original_df[selected_variable]
                cps_model = model_df["ngrdcol"] / model_df[selected_variable]

                abs_diff = np.abs( cps - cps_model )
                rms_error = np.sqrt( np.mean( abs_diff**2 ) )

                params_opt["rms_error"] = rms_error

                # Update table
                if gpu_fit:
                    gpu_batched_table_data.append({
                        col["id"]: f"{params_opt.get(col['id']):.3e}" if isinstance(params_opt.get(col["id"]), (int, float)) else params_opt.get(col["id"])
                        for col in gpu_batched_model_columns
                    })
                    gpu_batched_table_data[-1]["name"] = f"{filename}"
                else:
                    cpu_batched_table_data.append({
                        col["id"]: f"{params_opt.get(col['id']):.3e}" if isinstance(params_opt.get(col["id"]), (int, float)) else params_opt.get(col["id"])
                        for col in cpu_batched_model_columns
                    })
                    cpu_batched_table_data[-1]["name"] = f"{filename}"


                start_idx = end_idx

                
        if plot_mode == "cps":
            #fig = px.line(combined_df, x="Grid Boxes", y="Grid Boxes per Second", color="Name")
            #fig.update_layout( yaxis_title = "Throughput ( grid boxes / second)" )
            #fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Name")
            fig = px.line(
                combined_df, 
                x="ngrdcol", 
                y="Columns per Second", 
                color="Name", 
                symbol="Source"
            )
            fig.update_layout( 
                xaxis_title="Batch size (columns)", 
                yaxis_title = "Throughput (columns / second)"
            )
        else:
            #fig = px.line(combined_df, x="Grid Boxes", y=selected_variable, color="Name")
            fig = px.line(
                combined_df, 
                x="ngrdcol", 
                y=selected_variable, 
                color="Name", 
                symbol="Source"
            )
            fig.update_layout( yaxis_title = f"Runtime of {selected_variable} (seconds)" )

        for trace in fig.data:
            if trace.name.endswith("model"):
                trace.line.update(dash='dot')
                
        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            autosize=True,
            uirevision='constant',
            legend_title_text="Case"
        )
        fig.update_traces(marker=dict(size=8))
        fig = plot_with_enhancements(fig, "CPU Model: Throughput vs Batch Size")

        return fig.to_dict(), gpu_batched_table_data, cpu_batched_table_data


    # ======================================== App run ========================================
    app.run(host="0.0.0.0",debug=True,port=8051)



# ======================================== Main / helper functions ========================================
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

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Plot shared variables across multiple CSV files.")
    parser.add_argument("-dir", required=True, help="Directory to look for .csv files.")
    args = parser.parse_args()

    csv_files = glob(os.path.join(args.dir, "*.csv"))
    grouped_files = group_files_by_case(csv_files)

    all_variables, variable_sets = get_all_variables(
        [file for files in grouped_files.values() for file in files.values()]
    )

    launch_dash_app(args.dir,grouped_files, all_variables)
