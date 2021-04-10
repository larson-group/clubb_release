# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd

import numpy as np
import pdb

from analyze_sensitivity_matrix import analyzeSensMatrix
from test_analyzeSensMatrix import write_test_netcdf_files

# The parameters are tunable model parameters.
#    The order of paramsNames must match the order of filenames below.
paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

# Filenames of netcdf files containing fake testing data.
# Netcdf file containing metric and parameter values from the default simulation:
defaultNcFilename = 'default.nc'
#   sensNcFilenames lists filenames of sensitivity simulations.
#   Each file contains metrics and parameter values for a single simulation.
#   There should be one sensitivity simulation per each tunable parameter.
#   The 1st sens filename perturbs the first param in paramsNames, etc.
sensNcFilenames = np.array(['sens0.nc', 'sens1.nc'])

# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# Global mean, NOTE PRECT is in the unit of m/s
obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

# Write obsMetricValsDict to netcdf file in order to test whether SVD solution recovers original
#    observations.
obsNcFilename = 'obs.nc'

# The metrics are observed quantities that we want a tuned simulation to match.
#    The order of metricNames determines the order of rows in sensMatrix.
metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

# Column vector of (positive) weights.  A small value de-emphasizes
#   the corresponding metric in the fit.
metricsWeights = np.array([[1.0], [0.000001], [1.0]])

# This is the subset of paramsNames that vary from [0,1] (e.g., C5)
#    and hence will be transformed to [0,infinity] in order to make
#    the relationship between parameters and metrics more linear:
transformedParams = np.array([''])

# Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
clubbC8Vals = np.array([0.2, 0.7, 0.2])

# Create fake testing data with simple values
write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

# Calculate changes in parameter values needed to match metrics.
defaultBiasesCol, sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    dparamsSoln, paramsSoln, defaultBiasesApprox = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParams,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Plot the biases of the default simulation and the SVD approximation of that
biasesMatrix = np.dstack((defaultBiasesCol,defaultBiasesApprox)).squeeze()
df = pd.DataFrame(biasesMatrix,
                  index=metricsNames,
                  columns= ['defBias', 'defBiasApprox'])
biasesFig = px.line(df, x=df.index, y=df.columns,
              title = 'Biases of default simulation and approximations thereof.')
biasesFig.update_yaxes(title="Bias")
biasesFig.update_xaxes(title="Metric and region")

# Plot each column of normalized sensitivity matrix as a separate line.
# Each column tells us how all metrics vary with a single parameter.
df = pd.DataFrame(normlzdSensMatrix,
                  index=metricsNames,
                  columns=paramsNames)
normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized sensitivity matrix.<br>
                       Each column (line) shows how sensitive the metrics are<br>
                       to a change in a single parameter value.""" )
normlzdSensMatrixColsFig.update_yaxes(title="Normalized sensitivity, (1/metric) * dmetric/dparam")
normlzdSensMatrixColsFig.update_xaxes(title="Metric and region")

# Plot each row of normalized sensitivity matrix as a separate line.
# Each row tells us how a single metric varies with all parameters.
df = pd.DataFrame(np.transpose(normlzdSensMatrix),
                  index=paramsNames,
                  columns=metricsNames)
normlzdSensMatrixRowsFig = px.line(df, x=df.index, y=df.columns,
              title = """Rows of normalized sensitivity matrix.<br>
                       Each row (line) tells us how a single metric<br>
                       varies with all parameters.""" )
normlzdSensMatrixRowsFig.update_yaxes(title="Normalized sensitivity, (1/metric) * dmetric/dparam")
normlzdSensMatrixRowsFig.update_xaxes(title="Parameter")

# Plot the parameter values recommended by SVD.
paramsMatrix = np.dstack((dparamsSoln,paramsSoln)).squeeze()
df = pd.DataFrame(paramsMatrix,
                  index=paramsNames,
                  columns= ['dparamsSoln', 'paramsSoln'])
paramsFig = px.line(df, x=df.index, y=df.columns,
              title = 'Parameter values (and values - default) recommended by SVD.')
paramsFig.update_yaxes(title="Parameter value")
paramsFig.update_xaxes(title="Parameter Name")

sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),

        dcc.Graph( id='biasesFig', figure=biasesFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
        dcc.Graph( id='paramsFig', figure=paramsFig )

])

if __name__ == '__main__':
        sensMatrixDashboard.run_server(debug=True)
