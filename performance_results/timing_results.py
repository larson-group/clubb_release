#!/bin/python3

import argparse
import os
import pandas as pd
import re
from glob import glob
from dash import Dash, dcc, html, Input, Output, State, dash_table, ctx, ALL, no_update
import plotly.express as px
import numpy as np
from collections import defaultdict
from dash.dependencies import Input, Output, ALL
from dash.exceptions import PreventUpdate
import pyperclip
from scipy.optimize import minimize
import warnings
from itertools import chain
import plotly.graph_objects as go

# ======================================== GPU and CPU Labelling Code ========================================
gpu_names = ["A100", "V100", "H100", "MI250X"]
cpu_names = ["AMD7763", "Intel6430", "Intel6240"]

def is_gpu(name): return any(k in name for k in gpu_names)
def is_cpu(name): return any(k in name for k in cpu_names)

#base_colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#e377c2"]  # red, blue, green, purple, brown, pink
cpu_palette = [
    "blue",  # blue
    "red",  # vermillion
    "green",  # bluish green
    "purple",  # reddish purple
    "pink",  # sky blue
]

gpu_palette = [
    "darkorange",   # orange to contrast vermillion/red
    "seagreen",     # warmer green
    "mediumorchid", # lighter, warmer purple
    "deeppink",     # vivid pink
    "dodgerblue",   # brighter, slightly lighter blue
]

cpu_symbols = ["square", "circle", "cross", "triangle-up", "star"]
gpu_symbols = [s + "-open" for s in cpu_symbols]

def lighten(hex_color, factor=0.55):
    # Blend toward white; smaller factor -> lighter
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    r = int(r + (255 - r) * (1 - factor))
    g = int(g + (255 - g) * (1 - factor))
    b = int(b + (255 - b) * (1 - factor))
    return f"#{r:02x}{g:02x}{b:02x}"

def complementary(hex_color):
    """Return complementary color (RGB values inverted)."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    # Complementary = inverse in RGB space
    return f"#{(255 - r):02x}{(255 - g):02x}{(255 - b):02x}"
    
def complementary_and_darken(hex_color, darken_factor=0.85):
    """
    Return a complementary color (RGB inverted) and slightly darken it.
    """
    # Remove "#" and parse RGB
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)

    # Complementary color in RGB space
    r_c, g_c, b_c = 255 - r, 255 - g, 255 - b

    # Darken the complementary color
    r_d = int(r_c * darken_factor)
    g_d = int(g_c * darken_factor)
    b_d = int(b_c * darken_factor)

    return f"#{r_d:02x}{g_d:02x}{b_d:02x}"

def get_symbols_and_colors( combined_df ):

    series = list(combined_df["Configuration"].unique())
    cpu_series = [s for s in series if is_cpu(s)]
    gpu_series = [s for s in series if is_gpu(s)]
    others = [s for s in series if s not in cpu_series + gpu_series]

    color_map, symbol_map = {}, {}

    # CPUs: first 5 colors + solid symbols
    for i, name in enumerate(cpu_series):
        color_map[name]  = cpu_palette[i % len(cpu_palette)]
        symbol_map[name] = cpu_symbols[i % len(cpu_symbols)]

    # GPUs: last 5 colors + open symbols
    for i, name in enumerate(gpu_series):
        color_map[name]  = gpu_palette[i % len(gpu_palette)]
        symbol_map[name] = gpu_symbols[i % len(gpu_symbols)]

    # Any unclassified configs: cycle remaining CPU colors with a neutral symbol
    for j, name in enumerate(others):
        color_map[name]  = cpu_palette[(len(cpu_series) + j) % len(cpu_palette)]
        symbol_map[name] = "circle"

    # Legend ordering (CPUs first, then GPUs, then others)
    category_order = {"Configuration": cpu_series + gpu_series + others}

    return category_order, color_map, symbol_map

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

growing_container = {
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

cm_to_in = lambda cm: cm / 2.54

graph_config = {
    'toImageButtonOptions': {
        'format': 'png',
        'filename': 'custom_image',
        'height': 400,
        'width': 800,
        'scale': 3.125  # 300 DPI equivalent
    }
}


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
    N_v = np.floor( cols_per_core / flops_per_vop )

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20

    T_vcpu = ( T_r * N_s + T_v * N_v + b ) * cp_func( c, k, o, avg_array_size_MB )

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

def cpu_batched_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):
    params = params * param_scale
    _, rms_error = model_cpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )
    return rms_error
def model_cpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    # Model time
    m_ik, m_k, b, c, k, o = params

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20
 
    T_cpu =  ( m_ik * ngrdcol * N_vlevs  ) * cp_func( c, k, o, avg_array_size_MB ) + m_k * N_vlevs + b

    return T_cpu, rms_error( ngrdcol, runtime, T_cpu )


# ======================================== VCPU Batched ========================================

def vcpu_batched_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ):
    params = params * param_scale
    _, rms_error = model_vcpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func )
    return rms_error
def model_vcpu_batched_time(params, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func):

    # Model time
    m_vik, m_rik, m_k, b, c, k, o = params

    cols_per_core = ngrdcol / N_tasks
    flops_per_vop = N_vsize / N_prec
    
    # Number of columns using residual operations 
    #N_r = N_tasks * ( cols_per_core % flops_per_vop )
    N_r =  ngrdcol % (N_tasks * N_vsize / N_prec)

    avg_array_size_MB = ngrdcol * N_vlevs * (N_prec/8) / 2**20
 
    T_vcpu = (( m_vik * (ngrdcol - N_r ) + m_rik * N_r + m_k  ) * N_vlevs + b ) * cp_func( c, k, o, avg_array_size_MB )

    return T_vcpu, rms_error( ngrdcol, runtime, T_vcpu )


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

def gpu_batched_objective(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs ):
    _, rms_error = model_gpu_batched_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs )
    return rms_error
def model_gpu_batched_time(params, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs):

    m_ik, m_k, b = param_scale * params

    T_bgpu = m_ik * ngrdcol * N_vlevs + m_k * N_vlevs + b

    return T_bgpu, rms_error( ngrdcol, runtime, T_bgpu, mode="abs_cps" )


def model_throughputs(ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, model_version, bootstrap_samps=0, bootstrap_conf=90, selected_cp_funcs=None):

    fit_params_CI = None

    #============================== VCPU Model ==============================

    if model_version == "vcpu":

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]
        T_r_est = runtime[np.where(cols_per_core == 1)] - b_est
        T_v_est = runtime[np.where(cols_per_core == flops_per_vop)] - b_est

        param_scale = np.abs( [ T_r_est[0], T_v_est[0], b_est, 1.0, 1.0, 1.0 ] )

        initial_guess = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]

        bounds = [  (0,None),       # T_r_est
                    (0,None),       # T_v_est
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


        T_model, _ = model_vcpu_time(params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

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

        param_scale = np.abs( [ m_est, b_est, 1.0, 1.0, 1.0 ] )
        initial_guess = [ 1.0, 1.0, 1.0, 1.0, 1.0 ]

        bounds = [  (0,None),       # m_ik
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


        T_model, _ = model_cpu_time(params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best)

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

        T_model, _ = model_gpu_time(param_scale, params_opt, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs)

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

        T_model, _ = model_gpu_batched_time( params_opt, param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs )

        params_opt = param_scale * params_opt

        fit_params = {
            "m_ik": f"{params_opt[0]:.2e}", 
            "m_k": f"{params_opt[1]:.2e}",
            "b": f"{params_opt[2]:.2e}"
        }

        if bootstrap_samps > 0:

            # Bootstrap procedure
            params_bootstrap = []

            n_points = len(runtime)
            rng = np.random.default_rng()

            for _ in range(bootstrap_samps):
                indices = rng.integers(0, n_points, size=n_points)  # resample with replacement
                ngrdcol_bs = ngrdcol[indices]
                runtime_bs = runtime[indices]
                N_prec_bs = N_prec[indices]
                N_vlevs_bs = N_vlevs[indices]

                try:
                    result_bs = minimize(
                        gpu_batched_objective, 
                        initial_guess, 
                        args=(param_scale, ngrdcol_bs, runtime_bs, N_tasks, N_vsize, N_prec_bs, N_vlevs_bs ),
                        method='L-BFGS-B',
                        bounds=bounds
                    )
                    params_bootstrap.append(result_bs.x * param_scale)  # save rescaled params
                except Exception as e:
                    # Skip failed fits
                    print(f"Bootstrap fit failed: {e}")
                    continue
            
            params_bootstrap = np.array(params_bootstrap)

            # Now get 5th and 95th percentile (for 90% confidence interval)
            params_low = np.percentile(params_bootstrap, (100 - bootstrap_conf) / 2, axis=0)
            params_high = np.percentile(params_bootstrap, 100 - (100 - bootstrap_conf) / 2, axis=0)

            # Assemble fit_params dictionary
            fit_params_CI = {
                "m_ik": f"({params_low[0]:.2e}, {params_high[0]:.2e})",
                "m_k": f"({params_low[1]:.2e}, {params_high[1]:.2e})",
                "b": f"({params_low[2]:.2e}, {params_high[2]:.2e})"
            }

    #============================== CPU Batched Model ==============================

    elif model_version == "cpu_batched":

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

        T_model, _ = model_cpu_batched_time( params_opt, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best )

        # Assemble fit_params dictionary
        fit_params = {
            "m_ik": f"{params_opt[0]:.2e}",
            "m_k": f"{params_opt[1]:.2e}",
            "b": f"{params_opt[2]:.2e}",
            "c": f"{params_opt[3]:.2e}",
            "k": f"{params_opt[4]:.2e}",
            "o": f"{params_opt[5]:.2e}",
            "cp_func": cache_pen_best.__name__,
        }


        if bootstrap_samps > 0:

            # Bootstrap procedure
            params_bootstrap = []

            n_points = len(runtime)
            rng = np.random.default_rng()

            for _ in range(bootstrap_samps):
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

            # Now get 5th and 95th percentile (for 90% confidence interval)
            params_low = np.percentile(params_bootstrap, (100 - bootstrap_conf) / 2, axis=0)
            params_high = np.percentile(params_bootstrap, 100 - (100 - bootstrap_conf) / 2, axis=0)

            # Assemble fit_params dictionary
            fit_params_CI = {
                "m_ik": f"({params_low[0]:.2e}, {params_high[0]:.2e})",
                "m_k": f"({params_low[1]:.2e}, {params_high[1]:.2e})",
                "b": f"({params_low[2]:.2e}, {params_high[2]:.2e})",
                "c": f"({params_low[3]:.2e}, {params_high[3]:.2e})",
                "k": f"({params_low[4]:.2e}, {params_high[4]:.2e})",
                "o": f"({params_low[5]:.2e}, {params_high[5]:.2e})",
            }

    #============================== VCPU Batched Model ==============================
    elif model_version == "vcpu_batched":

        cols_per_core = ngrdcol / N_tasks
        flops_per_vop = N_vsize / N_prec

        b_est = runtime[0] - ( runtime[1] - runtime[0] ) / ( ngrdcol[1] - ngrdcol[0] ) * ngrdcol[0]
        m_k_est = b_est / N_vlevs[-1]
        first_vpoint = np.where(cols_per_core == flops_per_vop)[0][0]

        T_v_est = ( runtime[first_vpoint] - b_est ) / ( ngrdcol[first_vpoint] ) * ( 1 / N_vlevs[-1])
        T_r_est = ( runtime[0] - b_est ) / ( ngrdcol[1] ) * ( 1 / N_vlevs[-1])

        param_scale = np.abs( [ T_v_est, T_r_est, m_k_est, b_est, 1.0, 1.0, 1.0 ] )

        initial_guess = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ]

        bounds = [  (0,None),       # T_v
                    (0,None),       # T_r
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

            result = minimize(  vcpu_batched_objective, 
                                initial_guess, 
                                args=( param_scale, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cp_func ), 
                                method='CG' )
                                #method='L-BFGS-B',
                                #bounds = bounds )

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

        T_model, _ = model_vcpu_batched_time( params_opt, ngrdcol, runtime, N_tasks, N_vsize, N_prec, N_vlevs, cache_pen_best )

        # Assemble fit_params dictionary
        # fit_params = {
        #     "T_v": f"{params_opt[0]:.2e} - {param_scale[0]}",
        #     "T_r": f"{params_opt[1]:.2e} - {param_scale[1]}",
        #     "m_k": f"{params_opt[2]:.2e} - {param_scale[2]}",
        #     "b": f"{params_opt[3]:.2e} - {param_scale[3]}",
        #     "c": f"{params_opt[4]:.2e} - {param_scale[4]}",
        #     "k": f"{params_opt[5]:.2e} - {param_scale[5]}",
        #     "o": f"{params_opt[6]:.2e} - {param_scale[6]}",
        #     "cp_func": cache_pen_best.__name__,
        #     #"rms_error": rms_min
        # }

        fit_params = {
            "m_vik": f"{params_opt[0]:.2e}",
            "m_rik": f"{params_opt[1]:.2e}",
            "m_k": f"{params_opt[2]:.2e}",
            "b": f"{params_opt[3]:.2e}",
            "c": f"{params_opt[4]:.2e}",
            "k": f"{params_opt[5]:.2e}",
            "o": f"{params_opt[6]:.2e}",
            "cp_func": cache_pen_best.__name__,
            #"rms_error": rms_min
        }


        if bootstrap_samps > 0:

            # Bootstrap procedure
            params_bootstrap = []

            n_points = len(runtime)
            rng = np.random.default_rng()

            for _ in range(bootstrap_samps):
                indices = rng.integers(0, n_points, size=n_points)  # resample with replacement
                ngrdcol_bs = ngrdcol[indices]
                runtime_bs = runtime[indices]
                N_prec_bs = N_prec[indices]
                N_vlevs_bs = N_vlevs[indices]
                #N_tasks_bs = N_tasks if np.ndim(N_tasks) == 0 else N_tasks[indices]
                #N_vsize_bs = N_vsize if np.ndim(N_vsize) == 0 else N_vsize[indices]

                try:
                    result_bs = minimize(
                        vcpu_batched_objective, 
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

            # Now get 5th and 95th percentile (for 90% confidence interval)
            params_low = np.percentile(params_bootstrap, (100 - bootstrap_conf) / 2, axis=0)
            params_high = np.percentile(params_bootstrap, 100 - (100 - bootstrap_conf) / 2, axis=0)


            fit_params_CI = {
                "m_vik": f"({params_low[0]:.2e}, {params_high[0]:.2e})",
                "m_rik": f"({params_low[1]:.2e}, {params_high[1]:.2e})",
                "m_k": f"({params_low[2]:.2e}, {params_high[2]:.2e})",
                "b": f"({params_low[3]:.2e}, {params_high[3]:.2e})",
                "c": f"({params_low[4]:.2e}, {params_high[4]:.2e})",
                "k": f"({params_low[5]:.2e}, {params_high[5]:.2e})",
                "o": f"({params_low[6]:.2e}, {params_high[6]:.2e})",
            }

    else:

        print(f"ERROR: unknown fit mode {model_version} ")

    return T_model, fit_params, fit_params_CI


def launch_dash_app(dir_name, grouped_files, all_variables):

    data = {
        filename: pd.read_csv(filepath, comment="#")
        for group, files in grouped_files.items()
        for filename, filepath in files.items()
    }

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
            html.Details(
                [
                    html.Summary("Columns per Second Plot", style={"font-weight": "bold", "background-color": "#ddd", "padding": "10px", "cursor": "pointer", "font-size": "20px"}),
                    html.Div([

                        html.Div([
                            html.Label("Title: ", style={"margin-top": "5px"}),
                            dcc.Input(id="cps-plot-title", type="text", placeholder="Enter Title", value="Throughput vs Batch Size", style={"width": "80%"}),
                        ], style={"margin-bottom": "8px", "width": "600px"}),

                        html.Div([
                            html.Label("Legend Filter (Regex): ", style={"margin-top": "5px"}),
                            dcc.Input(id="cps-config-name-regex", type="text", placeholder="Name Filter", 
                                      value="_gptl|_derecho|_frontier|_casper|_arm", style={"width": "100%"}),
                        ], style={"margin-bottom": "8px", "width": "800px"}),

                        html.Div([
                            html.Label("Legend X"),
                            dcc.Slider(
                                id='cps-x-legend',
                                min=0,
                                max=1,
                                step=0.01,
                                value=0.05,  # default value
                                marks=None,         # hide all tick labels
                                tooltip={"placement": "bottom", "always_visible": True},
                                updatemode='drag'   # update the graph continuously while dragging
                            ),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}),

                        html.Div([
                            html.Label("Legend Y"),
                            dcc.Slider(
                                id='cps-y-legend',
                                min=0,
                                max=1,
                                step=0.01,
                                value=0.95,  # default value
                                marks=None,         # hide all tick labels
                                tooltip={"placement": "bottom", "always_visible": True},
                                updatemode='drag'   # update the graph continuously while dragging
                            ),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}),
                        
                        dcc.Graph(id="plot-columns-per-second", style=graph_style, config=graph_config),

                        html.Div([
                            html.Label("Base Batch Size: ", style={"margin-left": "5px"}),
                            dcc.Input(id="base-batch-size", type="number", value=0),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}), #"display": "none" makes this hidden

                        
                    ], style=plot_div_style )
                ],
                style=growing_container,
                open=True 
            ),

            # Time plot
            html.Details(
                [
                    html.Summary("Runtime Plot", style={"font-weight": "bold", "background-color": "#ddd", "padding": "10px", "cursor": "pointer", "font-size": "20px"}),
                    html.Div([

                        html.Div([
                            html.Label("Title: ", style={"margin-top": "5px"}),
                            dcc.Input(id="raw-plot-title", type="text", placeholder="Enter Title", value="Runtime vs Batch Size", style={"width": "80%"}),
                        ], style={"margin-bottom": "8px", "width": "600px"}),

                        html.Div([
                            html.Label("Legend Filter (Regex) ", style={"margin-top": "5px"}),
                            dcc.Input(id="raw-config-name-regex", type="text", placeholder="Name Filter", 
                                      value="_gptl|_derecho|_frontier|_casper|_arm", style={"width": "100%"}),
                        ], style={"margin-bottom": "8px", "width": "800px"}),

                        dcc.Graph(id="plot-raw", style=graph_style, config=graph_config),
                    
                        html.Summary("Display Modifications"),
                        dcc.RadioItems(
                            id="runtime-mods",
                            options=[
                                {"label": "none", "value": "none"},
                                {"label": "linear-fit", "value": "linfit"},
                                {"label": "scaling", "value": "scaling"}
                            ],
                            value="none",
                            labelStyle={"display": "inline-block"}  # 2 columns
                        ),
                        
                        html.Div([
                            dcc.RangeSlider(
                                id="range-slider",
                                min=1,
                                max=65536,
                                value=[128, 1024],  # Initial range selection
                                tooltip={"always_visible": True}
                            ),
                        ], style={"width": "60%", "padding": "20px"})
                        
                    ], style=plot_div_style ),
                ],
                style=growing_container,
                open=False
            ),

            # Fit plot
            html.Details(
                [
                    html.Summary("Individual Fit Plot", style={"font-weight": "bold", "background-color": "#ddd", "padding": "10px", "cursor": "pointer", "font-size": "20px"}),
                    html.Div([

                        # Plot
                        dcc.Graph(id="fit-plot", style=graph_style, config=graph_config), 

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
                ],
                style=growing_container,
                open=False
            ),

            # Batched fit plot
            html.Details(
                [
                    html.Summary("Batched Fit Plot", style={"font-weight": "bold", "background-color": "#ddd", "padding": "10px", "cursor": "pointer", "font-size": "20px"}),
                    html.Div([

                        html.Div([
                            html.Label("Title: ", style={"margin-top": "5px"}),
                            dcc.Input(id="fit-batched-plot-title", type="text", placeholder="Enter Title", value="Throughput vs Batch Size", style={"width": "80%"}),
                        ], style={"margin-bottom": "8px", "width": "600px"}),

                        html.Div([
                            html.Label("Legend Filter (Regex): ", style={"margin-top": "5px"}),
                            dcc.Input(id="fit-batched-config-name-regex", type="text", placeholder="Name Filter", 
                                      value="_gptl|_derecho|_frontier|_casper|_arm", style={"width": "100%"}),
                        ], style={"margin-bottom": "8px", "width": "800px"}),

                        html.Div([
                            html.Label("Legend X"),
                            dcc.Slider(
                                id='fit-batched-x-legend',
                                min=0,
                                max=1,
                                step=0.01,
                                value=0.05,  # default value
                                marks=None,         # hide all tick labels
                                tooltip={"placement": "bottom", "always_visible": True}
                            ),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}),

                        html.Div([
                            html.Label("Legend Y"),
                            dcc.Slider(
                                id='fit-batched-y-legend',
                                min=0,
                                max=1,
                                step=0.01,
                                value=0.95,  # default value
                                marks=None,         # hide all tick labels
                                tooltip={"placement": "bottom", "always_visible": True}
                            ),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}),

                        html.Div([
                            html.Label("Config Scale"),
                            dcc.Slider(
                                id='fit-batched-config-scale',
                                min=0,
                                max=2,
                                step=0.05,
                                value=1.2,  # default value
                                marks=None,         # hide all tick labels
                                tooltip={"placement": "bottom", "always_visible": True}
                            ),
                        ], style={"margin-bottom": "8px", "width": "200px", "display": "none"}),

                        # Plot
                        dcc.Graph(id="fit-plot-batched", style=graph_style, config=graph_config), 
                        
                        html.Label("Model Parameters"),
                        dash_table.DataTable(
                            id='batched-param-table',
                            #columns=gpu_batched_model_columns,
                            style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                            data=[],  # Initial empty table
                            editable=True,
                            row_deletable=False  # Allows deleting rows directly
                        ),

                        html.Label("RMS Error"),
                        dash_table.DataTable(
                            id='batched-rms-table',
                            #columns=gpu_batched_model_columns,
                            style_cell={'textAlign': 'center', "whiteSpace": "pre-line"},
                            data=[],  # Initial empty table
                            editable=True,
                            row_deletable=False  # Allows deleting rows directly
                        ),

                    ], style=fit_plot_div_style),
                ],
                style=growing_container,
                open=True
            ),

            # Custom function plot
            html.Details(
                [
                    html.Summary("Custom Function Plot", style={"font-weight": "bold", "background-color": "#ddd", "padding": "10px", "cursor": "pointer", "font-size": "20px"}),
                    
                    html.Div([

                        # Plot
                        dcc.Graph(id="plot-custom", style=graph_style, config=graph_config),

                        # Custom function input
                        dcc.Input(
                            id="custom-function",
                            type="text",
                            placeholder="Enter function (e.g., baseline + 5 * iter_total)",
                            debounce=True,
                            style={"width": "80%", "margin-top": "10px", "text-align": "center"}
                        ),

                    ], style=plot_div_style),
                ],
                style=growing_container,
                open=False
            ),

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
                        html.Summary(group.replace('_nz_', '_').replace('_', ' ')),
                        html.Div([
                            html.Button("Select All", id={"type": "select-all-button", "group": group}),
                            html.Button("Deselect All", id={"type": "deselect-all-button", "group": group}),
                            dcc.Checklist(
                                id={"type": "file-checkbox", "group": group},
                                options=[
                                    {"label": filename, "value": filename}
                                    for filename in filenames
                                ],
                                value=[
                                    filename
                                    for filename in list(filenames)[:2]
                                ] if i == 0 and filenames else []
                            )
                        ])
                    ],
                    open=(i == 0)
                )
                for i, (group, filenames) in enumerate(grouped_files.items())
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

            # Time / cps choice
            html.Div([
                dcc.RadioItems(
                    id="fit-plot-mode",
                    options=[
                        {"label": "Time", "value": "time"},
                        {"label": "Cols per Sec", "value": "cps"},
                        {"label": "Grid Boxes per Sec", "value": "gbps"}
                    ],
                    value="cps",
                    inline=True,
                ),
            ], style={"margin-bottom": "8px"}),

            # Fit function controls
            html.Div([
                
                html.Div([
                    html.Label("Fit Function Variables:", style={"font-weight": "bold"}),
                ], style={"margin-bottom": "8px"}),

                html.Div([

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
                        "Fit Individuals",
                        id="fit-function-button",
                        n_clicks=0,
                        style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer", "margin-bottom": "8px"}
                    ),

                    html.Div([
                        html.Label("Bootstrap Samples: ", style={"margin-top": "5px"}),
                        dcc.Input(id="bootstrap_samps", type="number", value=0, style={"width": "20%"}),
                    ], style={"margin-bottom": "8px"}),

                    html.Div([
                        html.Label("Bootstrap Confidence: ", style={"margin-top": "5px"}),
                        dcc.Input(id="bootstrap_conf", type="number", value=90, style={"width": "20%"}),
                    ], style={"margin-bottom": "8px"}),

                    # Batched fit button
                    html.Button(
                        "Fit CPU Batched",
                        id="fit-cpu-batched-function-button",
                        n_clicks=0,
                        style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
                    ),

                    html.Button(
                        "Fit VCPU Batched",
                        id="fit-vcpu-batched-function-button",
                        n_clicks=0,
                        style={"padding": "10px", "background-color": "#4CAF50", "color": "white", "border": "none", "cursor": "pointer"}
                    ),

                    html.Button(
                        "Fit GPU Batched",
                        id="fit-gpu-batched-function-button",
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


    def plot_with_enhancements(fig, title, scale_factor=1.2, x_legend=0.05, y_legend=0.95, config_scale=1):

        fig.update_traces(mode="lines+markers", selector=dict(mode="lines"))

        fig.update_traces(marker=dict(size=8))
        fig.update_traces(line=dict(width=2))
        
        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=18 * scale_factor)
            ),
            xaxis=dict(
                showgrid=True,
                gridcolor="lightgray",
                title_font=dict(size=14 * scale_factor),
                tickfont=dict(size=12 * scale_factor),
                showline=True,
                linewidth=1,
                linecolor="black",
                mirror=True
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor="lightgray",
                zeroline=True,
                range=[0, None],
                title_font=dict(size=14 * scale_factor),
                tickfont=dict(size=12 * scale_factor),
                showline=True,
                linewidth=1,
                linecolor="black",
                mirror=True
            ),
            legend=dict(
                x=0.01, y=0.99,
                xanchor="left", yanchor="top",
                font=dict(size=12 * config_scale),
                bgcolor="rgba(255,255,255,0.5)"
            ),
            # legend=dict(
            #     x=x_legend,          # 0 is left, 1 is right
            #     y=y_legend,          # 0 is bottom, 1 is top
            #     xanchor='left', # anchor the x position
            #     yanchor='top',   # anchor the y position
            #     bgcolor='rgba(255,255,255,0.5)',  # optional translucent background
            #     #bordercolor='black',
            #     #borderwidth=1,
            #     font=dict(size=12 * config_scale)
            # ),
            font=dict(size=12 * scale_factor),
            margin=dict(l=10, r=10, t=50, b=10),
            paper_bgcolor="white",
            plot_bgcolor="white",
        )

        return fig


    # ======================================== Time plot ========================================
    @app.callback(
    Output("plot-raw", "figure"),
        [
            Input({"type": "file-checkbox", "group": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("runtime-mods", "value"),
            Input("range-slider", "value"),
            Input("raw-plot-title", "value"),
            Input("raw-config-name-regex", "value")
        ]
    )
    def update_raw_plot(selected_files, selected_variable, xaxis_scale, yaxis_scale, runtime_mods, range_slider, title, config_regex ):

        # Flatten the list of lists into a single list of selected filenames
        selected_flat = list(chain.from_iterable(selected_files))

        if not selected_flat:
            raise PreventUpdate

        combined_df = pd.DataFrame()

        trendline_labels = []

        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:
                df = data[filename]
                if selected_variable in df.columns:  # Check before accessing

                    config_name = re.sub(config_regex, '', f"{filename}")

                    original_df = data[filename][["ngrdcol", selected_variable]].copy()
                    original_df = original_df[original_df["ngrdcol"] >= 32] 

                    original_df["Configuration"] = config_name
                    original_df["Source"] = "original"

                    if runtime_mods == "scaling":

                        if any(cpu in filename for cpu in cpu_names):

                            coeffs = np.polyfit(
                                original_df["ngrdcol"].where(original_df["ngrdcol"].between(*range_slider)).dropna(),
                                original_df[selected_variable].where(original_df["ngrdcol"].between(*range_slider)).dropna(),
                                deg=1
                            )

                            scale_curve = original_df.copy()

                            scale_curve[selected_variable] = scale_curve[selected_variable] / ( scale_curve["ngrdcol"] * coeffs[0] + coeffs[1] )
                            scale_curve["Source"] = "scale"

                            combined_df = pd.concat([combined_df, scale_curve])

                    else:

                        if runtime_mods == "linfit":

                            model_df = original_df.copy()

                            if any(cpu in filename for cpu in cpu_names):

                                # Use the slider to deremine CPU fit range
                                coeffs = np.polyfit(
                                    original_df["ngrdcol"].where(original_df["ngrdcol"].between(*range_slider)).dropna(),
                                    original_df[selected_variable].where(original_df["ngrdcol"].between(*range_slider)).dropna(),
                                    deg=1
                                )
                                label_position = 24
                            
                            else:

                                # Use all samples otherwise 
                                coeffs = np.polyfit(original_df["ngrdcol"], original_df[selected_variable], deg=1)
                                label_position = 11

                            model_df[selected_variable] = model_df["ngrdcol"] * coeffs[0] + coeffs[1]

                            model_df["Configuration"] = config_name
                            model_df["Source"] = "fit"

                            combined_df = pd.concat([combined_df, model_df])

                            # Create annotation text
                            m, b = coeffs
                            if any( cpu in config_name for cpu in cpu_names ):
                                equation_text = f"R\u0302<sub>CPU</sub>(N<sub>i</sub>) = {m:.3e}·N<sub>i</sub> + {b:.3e}"
                            else:
                                equation_text = f"R\u0302<sub>GPU</sub>(N<sub>i</sub>) = {m:.3e}·N<sub>i</sub> + {b:.3e}"

                            # Store one x,y point near end for label positioning
                            x_label = model_df["ngrdcol"][label_position]/1.3
                            #y_label = model_df[selected_variable][model_df["ngrdcol"] == x_label].values[0]
                            y_label = coeffs[1] + coeffs[0] * x_label

                            # Save label info to annotate later
                            trendline_labels.append({
                                "x": x_label,
                                "y": y_label,
                                "text": equation_text,
                                "config": config_name
                            })
                                
                        combined_df = pd.concat([combined_df, original_df])

        if combined_df.empty:
            raise PreventUpdate

            
        category_order, color_map_cfg, symbol_map_cfg = get_symbols_and_colors(combined_df)

        # 1) SeriesID: plain config for originals; suffix only for others
        combined_df["SeriesID"] = np.where(
            combined_df["Source"].eq("original"),
            combined_df["Configuration"],
            combined_df["Configuration"] + " • " + combined_df["Source"],
        )

        # 2) Create maps keyed by SeriesID
        symbol_map_series = {}
        dash_map_series   = {}
        color_map_series  = {}

        for cfg in combined_df["Configuration"].unique():
            base_color  = color_map_cfg.get(cfg, "black")
            base_symbol = symbol_map_cfg.get(cfg, "circle")

            # originals: key is just the config name (no ' • original')
            sid = cfg
            if (combined_df["SeriesID"] == sid).any():
                color_map_series[sid]  = base_color
                symbol_map_series[sid] = base_symbol
                dash_map_series[sid]   = "solid"

            # scale
            sid = f"{cfg} • scale"
            if (combined_df["SeriesID"] == sid).any():
                color_map_series[sid]  = base_color
                symbol_map_series[sid] = base_symbol
                dash_map_series[sid]   = "dash"

            # fit: black + circle + dotted
            sid = f"{cfg} • fit"
            if (combined_df["SeriesID"] == sid).any():
                color_map_series[sid]  = "black"
                symbol_map_series[sid] = "circle"
                dash_map_series[sid]   = "dot"
                # (if you later want no markers / no legend, handle post-plot)

        # 3) Order (originals → scales → fits) so originals draw underneath
        series_order = (
            [s for s in color_map_series if " • " not in s] +
            [s for s in color_map_series if " • scale" in s] +
            [s for s in color_map_series if " • fit"   in s]
        )

        # 4) Plot, legend title removed
        fig = px.line(
            combined_df,
            x="ngrdcol",
            y=selected_variable,
            color="SeriesID",
            symbol="SeriesID",
            line_dash="SeriesID",
            labels={"SeriesID": ""},  # no legend title
            category_orders={"SeriesID": series_order},
            color_discrete_map=color_map_series,
            symbol_map=symbol_map_series,
            line_dash_map=dash_map_series,
        )

        for trace in fig.data:
            if trace.name.endswith("fit"):
                trace.line.update(dash='dot')

        for label in trendline_labels:
            fig.add_annotation(
                x=label["x"],
                y=label["y"],
                text=label["text"],
                showarrow=True,
                arrowhead=1,
                ax=40,  # offset x for label placement
                ay=30 if any( cpu in label["config"] for cpu in cpu_names ) else -30,
                font=dict(size=16, color="black", family="Arial"),
                bgcolor="white",
                opacity=0.7
            )

        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            xaxis_title = "Batch Size (columns)", 
            yaxis_title = "Runtime (seconds)" if runtime_mods != "scaling" else "Cache Penalty Multiplier", 
            autosize=True,
            uirevision='constant'
        )

        fig = plot_with_enhancements(fig, title)

        for tr in fig.data:
            if "fit" in tr.name or "scale" in tr.name:
                tr.update(showlegend=False, mode="lines")

        return fig


    # ======================================== CPS plot ========================================
    @app.callback(
    Output("plot-columns-per-second", "figure"),
        [
            Input({"type": "file-checkbox", "group": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("cps-plot-title", "value"),
            Input("cps-config-name-regex", "value"),
            Input("base-batch-size", "value"),
            Input("cps-x-legend", "value"),
            Input("cps-y-legend", "value"),
            Input("fit-plot-mode", "value")
        ]
    )
    def update_columns_per_second_plot(selected_files, selected_variable, xaxis_scale, yaxis_scale, title, config_regex, base_batch_size, x_legend, y_legend, plot_mode):

        # Flatten the list of lists into a single list of selected filenames
        selected_flat = list(chain.from_iterable(selected_files))

        if not selected_flat:
            raise PreventUpdate

        combined_df = pd.DataFrame()

        ngrdcol_min = 0
        
        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:
                df = data[filename]
                if selected_variable in df.columns:  # Check before accessing

                    config_name = re.sub(config_regex, '', f"{filename}")

                    temp_df = data[filename].copy()

                    temp_df = temp_df[
                        (temp_df["ngrdcol"] >= max( 32, 2**int(base_batch_size) ) ) & 
                        (temp_df["ngrdcol"] % 2**int(base_batch_size) == 0)
                    ]

                    temp_df["ngrdcol"] = temp_df["ngrdcol"] / int( 2**int(base_batch_size) )

                    temp_df = temp_df[temp_df["ngrdcol"] >= temp_df["tasks"]]

                    if (ngrdcol_min is None) or (temp_df["ngrdcol"].min() > ngrdcol_min):
                        ngrdcol_min = temp_df["ngrdcol"].min()
                    

                    #print(temp_df["ngrdcol"])

                    temp_df["Configuration"] = config_name

                    
                    if plot_mode == "gbps":
                        temp_df["Throughput"] = temp_df["nz"] * temp_df["ngrdcol"] / temp_df[selected_variable]
                        temp_df["ngrdcol"] = temp_df["nz"] * temp_df["ngrdcol"]
                    else:
                        temp_df["Throughput"] = temp_df["ngrdcol"] / temp_df[selected_variable]

                    combined_df = pd.concat([combined_df, temp_df])

        if combined_df.empty:
            raise PreventUpdate
            
        combined_df = combined_df[combined_df["ngrdcol"] >= 32]

        category_order, color_map, symbol_map = get_symbols_and_colors( combined_df )

        fig = px.line(
            combined_df,
            x="ngrdcol",
            y="Throughput",
            color="Configuration",
            symbol="Configuration",
            labels={"Configuration": ""},
            category_orders=category_order,
            color_discrete_map=color_map,
            symbol_map=symbol_map,
        )

        # fig = px.line(combined_df, x="ngrdcol", y="Throughput", color="Configuration", symbol="Configuration", labels={"Configuration": ""}, 
        #               symbol_sequence=["circle", "square", "diamond", "cross", "x", "triangle-up", "star", "triangle-down", "bowtie", "hourglass"] )

        if plot_mode == "gbps":
            fig.update_layout(
                xaxis=dict(type=xaxis_scale),
                yaxis=dict(type=yaxis_scale),
                xaxis_title = "Batch size (grid boxes)",
                yaxis_title = "Throughput (grid boxes per second)",
                autosize=False,
                uirevision='constant'
            )
        else:
            fig.update_layout(
                xaxis=dict(type=xaxis_scale),
                yaxis=dict(type=yaxis_scale),
                xaxis_title = "Batch Size (columns)" if base_batch_size == 0 else f"Batch size (chunks of {2**base_batch_size} columns)", 
                yaxis_title = "Throughput (columns per second)" if base_batch_size == 0 else f"Throughput (column chunks per second)",
                autosize=False,
                uirevision='constant'
            )

        
        if (
            len(selected_flat) == 2 and
            any(any(gpu in f for gpu in gpu_names) for f in selected_flat) and
            any(any(cpu in f for cpu in cpu_names) for f in selected_flat)
        ):

            # Extract names
            name1, name2 = selected_flat
            config_name1 = re.sub(config_regex, '', f"{name1}")
            config_name2 = re.sub(config_regex, '', f"{name2}")

            if any(gpu_name in name1 for gpu_name in gpu_names ):
                gpu_name, cpu_name = config_name1, config_name2
            else:
                gpu_name, cpu_name = config_name2, config_name1

            cpu_df = combined_df[combined_df["Configuration"] == cpu_name]
            gpu_df = combined_df[combined_df["Configuration"] == gpu_name]

            # Find maxes
            cpu_Tmax = cpu_df["Throughput"].max()
            cpu_Tmax_i = cpu_df.loc[cpu_df["Throughput"].idxmax(), "ngrdcol"]
            gpu_Tmax = gpu_df["Throughput"].max()
            gpu_Tmax_i = gpu_df.loc[gpu_df["Throughput"].idxmax(), "ngrdcol"]


            # Upward vertical line from CPU peak to GPU Peak (PPR)
            fig.add_trace(go.Scatter(
                x=[min(cpu_Tmax_i, cpu_Tmax_i),max(cpu_Tmax_i, cpu_Tmax_i)],
                y=[min(cpu_Tmax, gpu_Tmax),max(cpu_Tmax, gpu_Tmax)],
                mode="lines+text",
                line=dict(color="black", dash="dot"),
                showlegend=False
            ))
            
            # Horizontal line atop PPR line
            fig.add_trace(go.Scatter(
                x=[cpu_Tmax_i/1.2, cpu_Tmax_i*1.2],
                y=[gpu_Tmax, gpu_Tmax],
                mode="lines+text",
                line=dict(color="black", dash="dot"),
                showlegend=False
            ))

            # Text for PPR
            fig.add_trace(go.Scatter(
                x=[cpu_Tmax_i*1.05],
                y=[(cpu_Tmax+gpu_Tmax)/2],
                mode="text",
                text=[f"PPR = {round_sigfigs(gpu_Tmax/cpu_Tmax, 3)}x"],
                textposition="middle right",
                textfont=dict(size=18, color="black"),
                showlegend=False
            ))

            cpu_final_val = cpu_df["Throughput"].values[-1]

            # Downward vertical line from GPU peak to CPU asymptote, ATR
            fig.add_trace(go.Scatter(
                x=[gpu_Tmax_i,gpu_Tmax_i],
                y=[gpu_Tmax,cpu_final_val],
                mode="lines+text",
                line=dict(color="black", dash="dot"),
                showlegend=False
            ))

            # Horizontal line atop ATR line
            fig.add_trace(go.Scatter(
                x=[gpu_Tmax_i/1.2, gpu_Tmax_i*1.2],
                y=[cpu_final_val, cpu_final_val],
                mode="lines+text",
                line=dict(color="black", dash="dot"),
                showlegend=False
            ))

            # Text for ATR
            fig.add_trace(go.Scatter(
                x=[gpu_Tmax_i/1.05],
                y=[(cpu_final_val+gpu_Tmax)/2],
                mode="text",
                text=[f"ATR = {round_sigfigs(gpu_Tmax/cpu_final_val, 3)}x"],
                textposition="middle left",
                textfont=dict(size=18, color="black"),
                showlegend=False
            ))



            # Find two points in dataset2 that straddle cpu_Tmax
            above = gpu_df[gpu_df["Throughput"] > cpu_Tmax]
            below = gpu_df[gpu_df["Throughput"] <= cpu_Tmax]

            # PRC line
            if not above.empty and not below.empty:
                # Take nearest points for interpolation
                y1 = below.iloc[-1]["Throughput"]
                x1 = below.iloc[-1]["ngrdcol"]
                y2 = above.iloc[0]["Throughput"]
                x2 = above.iloc[0]["ngrdcol"]

                # Linear interpolation to find crossing point
                x_cross = x1 + (cpu_Tmax - y1) * (x2 - x1) / (y2 - y1)

                # Do an inverse log scale transformation to find where this point is on the log plot
                x_cross_logged = 10.0**( np.log10(x1) + ( x_cross - x1 ) / ( x2 - x1 ) * ( np.log10(x2) - np.log10(x1) ) )

                # Horizontal line from CPU peak
                fig.add_trace(go.Scatter(
                    x=[cpu_Tmax_i, x_cross_logged],
                    y=[cpu_Tmax, cpu_Tmax],
                    mode="lines+text",
                    line=dict(color="black", dash="dot"),
                    showlegend=False
                ))

                # Downward vertical line from PRC
                fig.add_trace(go.Scatter(
                    x=[x_cross_logged, x_cross_logged],
                    y=[cpu_Tmax*0.03, cpu_Tmax],
                    mode="lines+text",
                    line=dict(color="black", dash="dot"),
                    showlegend=False
                ))

                # Text for PRC
                fig.add_trace(go.Scatter(
                    x=[x_cross],
                    y=[0],
                    mode="text",
                    text=[f"PRC = {int(round_sigfigs(x_cross, 3))}"],
                    textposition="bottom center",
                    textfont=dict(size=18, color="black"),
                    showlegend=False
                ))

            # Figure out EXACTLY where to place the DRC line
            for g_i in range(len(gpu_df)-1):

                gpu_i   = gpu_df.iloc[g_i]
                gpu_ip1 = gpu_df.iloc[g_i+1]

                for c_i in range(len(cpu_df)-1):

                    cpu_i   = cpu_df.iloc[c_i]
                    cpu_ip1 = cpu_df.iloc[c_i+1]

                    # Ensure the region is worth checking 
                    if ( 
                        max( gpu_i["ngrdcol"], cpu_i["ngrdcol"] ) < min( gpu_ip1["ngrdcol"], cpu_ip1["ngrdcol"] ) and
                        gpu_i["Throughput"] < cpu_i["Throughput"] and
                        gpu_ip1["Throughput"] > cpu_ip1["Throughput"]
                    ):

                        # We need to find the intersection between these points
                        xg1, yg1 = gpu_i["ngrdcol"], gpu_i["Throughput"]
                        xg2, yg2 = gpu_ip1["ngrdcol"], gpu_ip1["Throughput"]
                        xc1, yc1 = cpu_i["ngrdcol"], cpu_i["Throughput"]
                        xc2, yc2 = cpu_ip1["ngrdcol"], cpu_ip1["Throughput"]
                        
                        # Slopes and intercepts (assume no vertical lines and not parallel)
                        mg = (yg2 - yg1) / (xg2 - xg1)
                        bg = yg1 - mg * xg1

                        mc = (yc2 - yc1) / (xc2 - xc1)
                        bc = yc1 - mc * xc1

                        # Find x where mg*x + bg == mc*x + bc
                        x_int = (bc - bg) / (mg - mc)
                        y_int = mg * x_int + bg

                        # If the intersection of the interpolation is in the range we care about
                        if (min(xg1, xg2) <= x_int <= max(xg1, xg2) and
                            min(xc1, xc2) <= x_int <= max(xc1, xc2)):


                            def get_cross( xg1, xg2, xc1, xc2, yg1, yg2, yc1, yc2, do_log=False):

                                x_common = np.linspace(max(xc1, xg1), min(xc2, xg2), 5000)

                                def scaler( xt, x1, x2, do_log=False):
                                    if do_log:
                                        return ( np.log10(xt) - np.log10(x1) ) / ( np.log10(x2) - np.log10(x1) )
                                    else:
                                        return ( xt - x1 ) / ( x2 - x1 )

                                xg = xg1 + (xg2-xg1) * scaler( x_common, xg1, xg2, do_log )
                                xc = xc1 + (xc2-xc1) * scaler( x_common, xc1, xc2, do_log )

                                yg = yg1 + ( yg2 - yg1) * scaler( x_common, xg1, xg2, do_log )
                                yc = yc1 + ( yc2 - yc1) * scaler( x_common, xc1, xc2, do_log )

                                diff = yc - yg
                                cross_indices = np.where(np.diff(np.sign(diff)) != 0)[0]
                                
                                if len(cross_indices) > 0:
                                    i = cross_indices[0]

                                    # Linear interpolation of difference = 0
                                    d1 = yc[i] - yg[i]
                                    d2 = yc[i+1] - yg[i+1]
                                    alpha = -d1 / (d2 - d1)
                                    x_cross = x_common[i] + alpha * (x_common[i] - x_common[i+1])
                                    y_cross = yc[i] + alpha * (yc[i+1] - yc[i])

                                    return x_cross, y_cross

                                return None, None

                            # We need to transform these points to log world, THEN find the intersection
                            # this is because plotly connects points with linear line segments even when 
                            # plotting in log, so the x_cross we plot is not the actual one we want 
                            x_cross, y_cross = get_cross( xg1, xg2, xc1, xc2, yg1, yg2, yc1, yc2, do_log=False)
                            x_cross_adj, y_cross_adj = get_cross( xg1, xg2, xc1, xc2, yg1, yg2, yc1, yc2, do_log=True)

                            if x_cross_adj is None:
                                print("The x_cross transform broke the logic, we could move the transform up, but meh, defaulting to using linear x_cross")
                                x_cross_adj, y_cross_adj = x_cross, y_cross
                                print(f"x_cross = {x_cross}")
                
                            # Down ward vertical line from DRC
                            fig.add_trace(go.Scatter(
                                x=[x_cross_adj, x_cross_adj],
                                y=[y_cross_adj*0.05, y_cross_adj],
                                mode="lines+text",
                                line=dict(color="black", dash="dot"),
                                showlegend=False,
                                name="Intersection Marker"
                            ))

                            # Text for DRC
                            fig.add_trace(go.Scatter(
                                x=[x_cross],
                                y=[0],
                                mode="text",
                                text=[f"DRC = {int(round_sigfigs(x_cross, 3))}"],
                                textposition="bottom center",
                                textfont=dict(size=18, color="black"),
                                showlegend=False,
                                name="Intersection Marker"
                            ))


        return plot_with_enhancements(fig, title, x_legend=x_legend, y_legend=y_legend)


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
            Input({"type": "file-checkbox", "group": ALL}, "value"),
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
        combined_df = pd.DataFrame()


        for filename in selected_flat:
            if filename in data:
                temp_df = data[filename].copy()
                temp_df["Source"] = filename

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
            State({"type": "file-checkbox", "group": ALL}, "value"),
            State("variable-dropdown", "value")
        ]
    )
    def copy_csv_to_clipboard(n_clicks, selected_files, selected_variable):

        selected_flat = list(chain.from_iterable(selected_files))

        # The order in 'header' will match the order of these files
        header = ["ngrdcol"] + selected_flat

        data_rows = defaultdict(dict)

        # Populate data_rows
        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:
                df = data[filename][["ngrdcol", selected_variable]].copy()
                for ngrdcol_val in df["ngrdcol"].unique():
                    if ngrdcol_val not in data_rows:
                        data_rows[ngrdcol_val] = {}
                    for idx, row in df[df["ngrdcol"] == ngrdcol_val].iterrows():
                        data_rows[ngrdcol_val][filename] = row.get(selected_variable, "")

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
            State({"type": "file-checkbox", "group": ALL}, "value"),
            Input("variable-dropdown", "value"),
            Input("x-axis-scale", "value"),
            Input("y-axis-scale", "value"),
            Input("fit-plot-mode", "value"),
            State({"type": "cp_func_checkbox"}, "value"),
        ],
    )
    def update_fit_plot(n_clicks, N_vsize, selected_files, selected_variable, xaxis_scale, yaxis_scale, plot_mode, selected_cp_funcs):
        
        selected_flat = list(chain.from_iterable(selected_files))

        if n_clicks == 0 or not selected_files:
            raise PreventUpdate

        combined_df = pd.DataFrame()

        vcpu_table_data = []
        cpu_table_data = []
        gpu_table_data = []

        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:

                original_df = data[filename][["ngrdcol", selected_variable]].copy()
                original_df["Name"] = f"{filename}"
                original_df["Source"] = "original"
                original_df["Columns per Second"] = original_df["ngrdcol"] / original_df[selected_variable]

                combined_df = pd.concat([combined_df, original_df])
                
                # Get the number of cores used 
                N_tasks = data[filename]["tasks"][0]
                N_vlevs = data[filename]["nz"][0]

                # Apply model_throughput to the selected variable and add it as "_dup"
                vcpu_df = original_df.copy()
                cpu_df = original_df.copy()
                gpu_model_df = original_df.copy()

                if "_sp_" in filename:
                    N_prec = 32
                else:
                    N_prec = 64
                
                if any(gpu_name in filename for gpu_name in gpu_names ):

                    print(f"modeling {filename} with gpu")
                    T_gpu, gpu_params, _ = model_throughputs(  gpu_model_df["ngrdcol"].values, 
                                                            gpu_model_df[selected_variable].values, 
                                                            N_tasks, 
                                                            N_vsize, 
                                                            N_prec, 
                                                            N_vlevs, 
                                                            "gpu")

                    gpu_model_df[selected_variable] = T_gpu
                    gpu_model_df["Columns per Second"] = gpu_model_df["ngrdcol"] / gpu_model_df[selected_variable]
                    gpu_model_df["Name"] = filename
                    gpu_model_df["Source"] = "gmodel"

                    gpu_table_data.append({
                        col["id"]: f"{gpu_params.get(col['id']):.3e}" if isinstance(gpu_params.get(col["id"]), (int, float)) else gpu_params.get(col["id"])
                        for col in gpu_model_columns
                    })
                    gpu_table_data[-1]["name"] = filename

                    combined_df = pd.concat([combined_df, gpu_model_df])

                else:

                    print(f"modeling {filename} with vcpu")
                    T_vcpu, vcpu_params, _ = model_throughputs( vcpu_df["ngrdcol"].values, 
                                                             vcpu_df[selected_variable].values, 
                                                             N_tasks, 
                                                             N_vsize, 
                                                             N_prec, 
                                                             N_vlevs, 
                                                             "vcpu",
                                                             selected_cp_funcs=selected_cp_funcs )

                    vcpu_df[selected_variable] = T_vcpu
                    vcpu_df["Columns per Second"] = vcpu_df["ngrdcol"] / vcpu_df[selected_variable]
                    vcpu_df["Name"] = filename
                    vcpu_df["Source"] = "vmodel"

                    vcpu_table_data.append({
                        col["id"]: f"{vcpu_params.get(col['id']):.3e}" if isinstance(vcpu_params.get(col["id"]), (int, float)) else vcpu_params.get(col["id"])
                        for col in vcpu_model_columns
                    })
                    vcpu_table_data[-1]["name"] = filename


                    print(f"modeling {filename} with cpu")
                    T_cpu, cpu_params, _ = model_throughputs( cpu_df["ngrdcol"].values, 
                                                           cpu_df[selected_variable].values, 
                                                           N_tasks, 
                                                           N_vsize, 
                                                           N_prec, 
                                                           N_vlevs, 
                                                           "cpu",
                                                           selected_cp_funcs=selected_cp_funcs)
                    
                    cpu_df[selected_variable] = T_cpu
                    cpu_df["Columns per Second"] = cpu_df["ngrdcol"] / cpu_df[selected_variable]
                    cpu_df["Name"] = filename
                    cpu_df["Source"] = "cmodel"

                    cpu_table_data.append({
                        col["id"]: f"{cpu_params.get(col['id']):.3e}" if isinstance(cpu_params.get(col["id"]), (int, float)) else cpu_params.get(col["id"])
                        for col in cpu_model_columns
                    })
                    cpu_table_data[-1]["name"] = filename


                    combined_df = pd.concat([combined_df, vcpu_df, cpu_df])
                
        if plot_mode == "cps":
            fig = px.line(combined_df, x="ngrdcol", y="Columns per Second", color="Name", symbol="Source")
            fig.update_layout( yaxis_title = "Throughput (columns / second)" )
        else:
            fig = px.line(combined_df, x="ngrdcol", y=selected_variable, color="Name", symbol="Source")
            fig.update_layout( yaxis_title = f"Runtime of {selected_variable} (seconds)" )


        for trace in fig.data:
            if trace.name.endswith("model"):
                trace.line.update(dash='dot')

        fig.update_layout(
            xaxis=dict(type=xaxis_scale),
            yaxis=dict(type=yaxis_scale),
            xaxis_title = "Batch Size (columns)", 
            autosize=True,
            uirevision='constant'
        )
        fig = plot_with_enhancements(fig, f"{selected_variable} vs. Number of Grid Columns (Raw Values)")

        return fig.to_dict(), vcpu_table_data, cpu_table_data, gpu_table_data


    # ======================================== Batched fit plot ========================================
    @app.callback(
        [
            Output("fit-plot-batched", "figure"),
            Output('batched-param-table', 'data'),
            Output('batched-rms-table', 'data'),
        ],
        [
            Input("fit-cpu-batched-function-button", "n_clicks"),
            Input("fit-vcpu-batched-function-button", "n_clicks"),
            Input("fit-gpu-batched-function-button", "n_clicks"),
            State("N_vsize", "value"),
            State("bootstrap_samps", "value"),
            State("bootstrap_conf", "value"),
            State({"type": "file-checkbox", "group": ALL}, "value"),
            State("variable-dropdown", "value"),
            State("x-axis-scale", "value"),
            State("y-axis-scale", "value"),
            State("fit-plot-mode", "value"),
            State({"type": "cp_func_checkbox"}, "value"),
            State("fit-batched-plot-title", "value"),
            State("fit-batched-config-name-regex", "value"),
            State("fit-batched-x-legend", "value"),
            State("fit-batched-y-legend", "value"),
            State("fit-batched-config-scale", "value")
        ],
        prevent_initial_call=True
    )
    def update_fit_plot_batched(cpu_clicks, vcpu_clicks, gpu_clicks, N_vsize, bootstrap_samps, bootstrap_conf, 
                                selected_files, selected_variable, xaxis_scale, yaxis_scale, 
                                plot_mode, selected_cp_funcs, title, config_regex, x_legend, y_legend, config_scale):

        cpu_fit, vcpu_fit, gpu_fit = False, False, False

        if ctx.triggered_id ==  'fit-cpu-batched-function-button':
            cpu_fit         = True 
            model_version   = "cpu_batched"
        elif ctx.triggered_id == 'fit-vcpu-batched-function-button':
            vcpu_fit         = True 
            model_version   = "vcpu_batched"
        elif ctx.triggered_id == 'fit-gpu-batched-function-button':
            gpu_fit         = True 
            model_version   = "gpu_batched"
            
        selected_flat = list(chain.from_iterable(selected_files))

        combined_df = pd.DataFrame()

        N_tasks = np.empty(0)
        N_vlevs = np.empty(0)
        runtimes = np.empty(0)
        ngrdcols = np.empty(0)
        N_prec = np.empty(0)
        filenames = ""


        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:
                
                original_df = data[filename].copy()
                original_df = original_df[original_df["ngrdcol"] > 32]

                # Get the number of cores used 
                N_tasks  = np.concatenate( [N_tasks, np.array( original_df["tasks"] )] )
                N_vlevs  = np.concatenate( [N_vlevs, np.array( original_df["nz"] )] )
                runtimes = np.concatenate( [runtimes, np.array( original_df[selected_variable] )] )
                ngrdcols = np.concatenate( [ngrdcols, np.array( original_df["ngrdcol"] )] )

                if "_sp_" in filename:
                    N_prec = np.concatenate( [N_prec, np.array([32]*len(original_df["nz"])) ])
                else:
                    N_prec = np.concatenate( [N_prec, np.array([64]*len(original_df["nz"])) ])
                

        T_model, params_opt, params_opt_CI = model_throughputs( ngrdcols, 
                                                 runtimes, 
                                                 N_tasks, 
                                                 N_vsize, 
                                                 N_prec, 
                                                 N_vlevs, 
                                                 model_version,
                                                 bootstrap_samps, 
                                                 bootstrap_conf,
                                                 selected_cp_funcs )

        batched_param_table = []

        # First row: params_opt
        row_data = {"": "Optimal Parameters", **{key: value for key, value in params_opt.items()}}
        batched_param_table.append(row_data)

        # Second row: params_opt_CI
        if params_opt_CI is not None:
            ci_row_data = {"": f"{bootstrap_conf}% Confidence Interval", **{key: value for key, value in params_opt_CI.items()}}
            batched_param_table.append(ci_row_data)
        
        batched_rms_table = []

        start_idx = 0
        
        for filename in selected_flat:
            if filename in data and selected_variable in data[filename].columns:

                config_name = re.sub(config_regex, '', f"{filename}")

                original_df = data[filename][["ngrdcol", selected_variable]].copy()
                original_df = original_df[original_df["ngrdcol"] > 32]

                original_df["Name"] = config_name
                original_df["Source"] = "original"
                original_df["Columns per Second"] = original_df["ngrdcol"] / original_df[selected_variable]

                end_idx = start_idx + len(original_df["ngrdcol"])

                model_df = original_df.copy()
                model_df["Name"] = config_name
                model_df["Source"] = "model"
                model_df[selected_variable] = T_model[start_idx:end_idx]
                model_df["Columns per Second"] = model_df["ngrdcol"] / model_df[selected_variable]
                
                combined_df = pd.concat([combined_df, original_df, model_df])
                start_idx = end_idx

                # GPU uses absolute error 
                rms_error_model = rms_error( model_df["ngrdcol"], 
                                             original_df[selected_variable], 
                                             model_df[selected_variable], 
                                             mode="abs_cps" if gpu_fit else "abs_percent_cps" )

                # Add name field first
                row_data = {"Configuration": f"{filename}" }
                row_data = {**row_data, "RMS Error": f"{rms_error_model}"}

                # Add row to table
                batched_rms_table.append(row_data)


                
        if plot_mode == "cps":
            fig = px.line(
                combined_df, 
                x="ngrdcol", 
                y="Columns per Second", 
                color="Name", 
                symbol="Name",
                line_dash="Source"
            )
            fig.update_layout( 
                xaxis_title="Batch Size (columns)", 
                yaxis_title = "Throughput (columns / second)"
            )
        else:
            fig = px.line(
                combined_df, 
                x="ngrdcol", 
                y=selected_variable, 
                color="Name", 
                symbol="Name",
                line_dash="Source"
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
            legend_title_text="Configuration"
        )

        #fig.update_traces(marker=dict(size=8))

        fig = plot_with_enhancements(fig, title, x_legend=x_legend, y_legend=y_legend, config_scale=config_scale)

        return fig.to_dict(), batched_param_table, batched_rms_table


    @app.callback(
        Output({"type": "file-checkbox", "group": ALL}, "value"),
        Input({"type": "select-all-button", "group": ALL}, "n_clicks"),
        Input({"type": "deselect-all-button", "group": ALL}, "n_clicks"),
        State({"type": "file-checkbox", "group": ALL}, "options"),
        prevent_initial_call=True
    )
    def select_deselect_all(select_clicks, deselect_clicks, options_lists):
        triggered = ctx.triggered_id
        if triggered is None:
            return [no_update] * len(options_lists)

        group_triggered = triggered["group"]
        type_triggered = triggered["type"]

        new_values = []
        for i, options in enumerate(options_lists):
            group_this = ctx.inputs_list[0][i]["id"]["group"]
            if group_this != group_triggered:
                new_values.append(no_update)
            elif type_triggered == "select-all-button":
                new_values.append([opt["value"] for opt in options])
            elif type_triggered == "deselect-all-button":
                new_values.append([])
            else:
                new_values.append(no_update)

        return new_values


    # ======================================== App run ========================================
    app.run(host="0.0.0.0",debug=True,port=8051)



# ======================================== Main / helper functions ========================================
def round_sigfigs(x, sigfigs=3):
    if x == 0:
        return 0
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (sigfigs - 1))

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

# def group_files_by_config(csv_files):
#     configs = defaultdict(dict)
    
#     for file in csv_files:
#         basename = os.path.basename(file)
#         # Match everything before ".csv"
#         match = re.match(r"^(.*)\.csv$", basename)
#         if match:
#             filename = match.group(1)

#             # Remove the "_<number>nz_" part to form the group key
#             group_key = re.sub(r'_\d+nz_', '_nz_', filename)

#             configs[group_key][filename] = file

#     # Sort filenames naturally within each config group
#     sorted_configs = {}
#     for config, files in configs.items():
#         sorted_configs[config] = dict(sorted(files.items(), key=lambda item: natural_key(item[0])))

#     return sorted_configs

def group_files_by_config(csv_files):
    configs = defaultdict(dict)
    
    for file in csv_files:
        basename = os.path.basename(file)
        # Match everything before ".csv"
        match = re.match(r"^(.*)\.csv$", basename)
        if match:
            filename = match.group(1)

            # Extract prefix before first underscore (e.g., "A100" from "A100_nvhpc_xyz")
            group_key_match = re.match(r"^([^_]+)", filename)
            if group_key_match:
                group_key = group_key_match.group(1)
                configs[group_key][filename] = file

    # Sort filenames naturally within each config group
    sorted_configs = {}
    for config, files in configs.items():
        sorted_configs[config] = dict(sorted(files.items(), key=lambda item: natural_key(item[0])))

    return sorted_configs

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Plot shared variables across multiple CSV files.")
    parser.add_argument(
        "directory",
        nargs="?",                 # makes it optional
        default="timing_results",  # default if not provided
        help="Directory to look for .csv files (default: timing_results)."
    )
    args = parser.parse_args()

    csv_files = glob(os.path.join(args.directory, "*.csv"))
    grouped_files = group_files_by_config(csv_files)

    all_variables, variable_sets = get_all_variables(
        [file for files in grouped_files.values() for file in files.values()]
    )

    launch_dash_app(args.directory, grouped_files, all_variables)
