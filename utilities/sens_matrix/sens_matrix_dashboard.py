# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

import numpy as np
import pdb

from analyze_sensitivity_matrix import analyzeSensMatrix
from test_analyzeSensMatrix import write_test_netcdf_files

# The parameters are tunable model parameters.
#    The order of paramsNames must match the order of filenames below.
#paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

# Filenames of netcdf files containing fake testing data.
# Netcdf file containing metric and parameter values from the default simulation:
#defaultNcFilename = 'default.nc'
#   sensNcFilenames lists filenames of sensitivity simulations.
#   Each file contains metrics and parameter values for a single simulation.
#   There should be one sensitivity simulation per each tunable parameter.
#   The 1st sens filename perturbs the first param in paramsNames, etc.
#sensNcFilenames = np.array(['sens0.nc', 'sens1.nc'])

# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# Global mean, NOTE PRECT is in the unit of m/s
#obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

# Write obsMetricValsDict to netcdf file in order to test whether SVD solution recovers original
#    observations.
#obsNcFilename = 'obs.nc'

# The metrics are observed quantities that we want a tuned simulation to match.
#    The order of metricNames determines the order of rows in sensMatrix.
#metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

# Column vector of (positive) weights.  A small value de-emphasizes
#   the corresponding metric in the fit.
#metricsWeights = np.array([[1.0], [0.000001], [1.0]])

# This is the subset of paramsNames that vary from [0,1] (e.g., C5)
#    and hence will be transformed to [0,infinity] in order to make
#    the relationship between parameters and metrics more linear:
#transformedParamsNames = np.array([''])

# Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
#clubbC8Vals = np.array([0.2, 0.7, 0.2])

# Create fake testing data with simple values
#write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

# Metrics are observed quantities that we want a tuned simulation to match.
#    The order of metricNames determines the order of rows in sensMatrix.
# Column vector of (positive) weights.  A small value de-emphasizes
#   the corresponding metric in the fit.
metricsNamesAndWeights = [ \
                        ['SWCF_GLB', 2.], \
                        ['SWCF_DYCOMS', 1.], \
                        ['SWCF_HAWAII', 1.], \
                        ['SWCF_VOCAL', 1.], \
                        ['SWCF_LBA', 1.], \
                        ['SWCF_WP', 1.], \
                        ['SWCF_EP', 1.], \
                        ['SWCF_NP', 1.], \
                        ['SWCF_SP', 1.], \
                        ['SWCF_PA', 1.], \
                        ['SWCF_CAF', 1.], \
                        ['LWCF_GLB', 2.], \
##                        ['LWCF_DYCOMS', 0.1], \
##                        ['LWCF_HAWAII', 0.1], \
##                        ['LWCF_VOCAL', 0.1], \
                        ['LWCF_LBA', 1.], \
                        ['LWCF_WP', 1.], \
                        ['LWCF_EP', 1.], \
                        ['LWCF_NP', 1.], \
                        ['LWCF_SP', 1.], \
                        ['LWCF_PA', 1.], \
                        ['LWCF_CAF', 1.], \
                        ['PRECT_GLB', 2.], \
##                        ['PRECT_DYCOMS', 0.1], \
##                        ['PRECT_HAWAII', 0.1], \
##                        ['PRECT_VOCAL', 0.1], \
                        ['PRECT_LBA', 1.], \
                        ['PRECT_WP', 1.], \
                        ['PRECT_EP', 1.], \
                        ['PRECT_NP', 1.], \
                        ['PRECT_SP', 1.], \
                        ['PRECT_PA', 1.], \
                        ['PRECT_CAF', 1.] \
                         ]


dfMetricsNamesAndWeights =  \
        pd.DataFrame( metricsNamesAndWeights, columns = ['metricsNames', 'metricsWeights'] )
metricsNames = dfMetricsNamesAndWeights[['metricsNames']].to_numpy().astype(str)[:,0]
metricsWeights = dfMetricsNamesAndWeights[['metricsWeights']].to_numpy()


# Parameters are tunable model parameters.
# Each parameter is associated with a single sensitivity simulation in which that parameter is perturbed,
#    and the output from each sensitivity simulation is expected to be stored in its own netcdf file.
#    Each netcdf file contains metrics and parameter values for a single simulation.
# paramsNames designates the parameter whose value is changed in the corresponding sensitivity simulation.
# The float listed below is a factor that is used below for scaling plots.
paramsNamesScalesAndFilenames = [ \
#                  ['clubb_c8', 1.0, 'devel/anvil.devel.C82.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.0, 'devel/anvil.devel.n2thres3.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_n2', 1.0, 'devel/anvil.devel.n21.ne30_ne30_Regional.nc'], \
#                  ['prc_exp', 1.0, 'devel/anvil.devel.prc2.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_wpxp_ri', 1.0, 'devel/anvil.devel.ri4.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_n2_clear_wp3', 1.0, 'devel/anvil.devel.wp35.ne30_ne30_Regional.nc'], \
                  ['clubb_c8', 1.0, '20210505/anvil.devel.0503_c8p4.ne30_ne30_Regional.nc'], \
                  ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210505/anvil.devel.0503_thresh2.2.ne30_ne30_Regional.nc'], \
                  ['clubb_c_invrs_tau_n2', 1.0, '20210505/anvil.devel.0503_n2p4.ne30_ne30_Regional.nc'], \
                  ['clubb_c_k10', 1.0, '20210505/anvil.devel.0503_ck10p3.ne30_ne30_Regional.nc'], \
                  ['vqit', 1.0, '20210505/anvil.devel.0503_vqitp4.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_wpxp_ri', 1.0, '20210505/anvil.devel.0503_ri5.ne30_ne30_Regional.nc'], \
#                  ['clubb_gamma_coef', 1.0, '20210505/anvil.devel.0503_gap2.ne30_ne30_Regional.nc'], \
#                  ['max_total_ni', 1.0, '20210505/anvil.devel.0503_ni700.ne30_ne30_Regional.nc'], \
#                  ['cldfrc_rhminl', 1.0, '20210505/anvil.devel.0503_rhp9.ne30_ne30_Regional.nc'] \
                        ]

dfparamsNamesScalesAndFilenames =  \
        pd.DataFrame( paramsNamesScalesAndFilenames, columns = ['paramsNames', 'paramsScales', 'sensNcFilenames'] )
paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:,0]
# Extract scaling factors of parameter values.
# The scaling allows us to avoid plotting very large or small values.
paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:,0]
sensNcFilenames = dfparamsNamesScalesAndFilenames[['sensNcFilenames']].to_numpy().astype(str)[:,0]

# This the subset of paramsNames that vary from [0,1] (e.g., C5)
#    and hence will be transformed to [0,infinity] in order to make
#    the relationship between parameters and metrics more linear:
transformedParamsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2', 'clubb_c_invrs_tau_n2_clear_wp3'])
#transformedParamsNames = np.array([''])

# Netcdf file containing metric and parameter values from the default simulation
defaultNcFilename = \
    '20210505/anvil.devel.NGD_thres2p3_vqitp5.ne30_ne30_Regional.nc'#'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux.ne30_ne30_GLBmean.nc'
#    'devel/anvil.devel.base.ne30_ne30_Regional.nc'

# Metrics from simulation that use the SVD-recommended parameter values
# Here, we use default simulation just as a placeholder.
linSolnNcFilename = \
   '20210505/anvil.devel.NGD_thres2p3_vqitp5.ne30_ne30_Regional.nc' #'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux.ne30_ne30_GLBmean.nc'

# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# Global mean, NOTE PRECT is in the unit of m/s
obsMetricValsDict = { \
   'LWCF_GLB': 28.008, 'PRECT_GLB': 0.000000031134259, 'SWCF_GLB': -45.81, 'TMQ_GLB': 24.423, \
   'LWCF_DYCOMS': 17.39229, 'PRECT_DYCOMS': 5.540433451401774e-09, 'SWCF_DYCOMS': -63.05767, 'TMQ_DYCOMS': 19.03481,\
   'LWCF_LBA': 54.9896, 'PRECT_LBA': 8.75291e-08, 'SWCF_LBA': -67.0833, 'TMQ_LBA': 48.1414,\
   'LWCF_HAWAII': 21.111, 'PRECT_HAWAII': 1.47959e-08, 'SWCF_HAWAII': -31.6793, 'TMQ_HAWAII': 30.7126,\
   'LWCF_WP': 56.3396, 'PRECT_WP': 7.65602e-08, 'SWCF_WP': -64.7297, 'TMQ_WP':50.7996,\
   'LWCF_EP': 32.724, 'PRECT_EP': 5.93702e-08, 'SWCF_EP': -54.0756, 'TMQ_EP':45.7974,\
   'LWCF_NP': 30.3083, 'PRECT_NP': 3.657e-08, 'SWCF_NP': -73.3256, 'TMQ_NP':12.9728,\
   'LWCF_SP': 30.6713, 'PRECT_SP': 3.11686e-08, 'SWCF_SP': -65.3633, 'TMQ_SP':12.5768,\
   'LWCF_PA': 28.5308, 'PRECT_PA': 3.22142e-08, 'SWCF_PA': -73.4048, 'TMQ_PA':50.0178,\
   'LWCF_CAF': 57.0782, 'PRECT_CAF': 5.71253e-08, 'SWCF_CAF': -60.8234, 'TMQ_CAF':41.1117,\
   'LWCF_VOCAL': 16.2196, 'PRECT_VOCAL': 1.78555e-09, 'SWCF_VOCAL': -77.2623, 'TMQ_VOCAL':17.5992  }

# Calculate changes in parameter values needed to match metrics.
defaultMetricValsCol, defaultBiasesCol, \
defaultBiasesApprox, defaultBiasesApproxPC, \
defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
defaultParamValsOrigRow, dparamsSoln, \
paramsSoln, paramsLowVals, paramsHiVals, \
paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
         analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Plot the biases of the default simulation and the SVD approximation of that
biasesMatrix = np.dstack((defaultBiasesCol,defaultBiasesApprox, defaultBiasesApproxPC)).squeeze()
fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(defaultMetricValsCol))) @ biasesMatrix
df = pd.DataFrame(fracBiasesMatrix,
                  index=metricsNames,
                  columns= ['fracDefBias', 'fracDefBiasApprox', 'defaultBiasesApproxPC'])
biasesFig = px.line(df, x=df.index, y=df.columns,
              title = 'Fractional biases of default simulation and approximations thereof.')
biasesFig.update_yaxes(title="-(Sim-Obs) / abs(default metric value)")
biasesFig.update_xaxes(title="Metric and region")
biasesFig.update_layout(hovermode="x")

# Compute errors in the final parameter solution as a frac of the original errors
# estBiasesOrig = default metric value + correction - observed metric value
defaultBiasesColSqd = defaultBiasesCol**2
#defaultMetricValsColSqd = defaultMetricValsCol**2
estBiasesSqd = (defaultBiasesApprox - defaultBiasesCol)**2
fracError = estBiasesSqd / defaultBiasesColSqd
#fracError = estBiasesSqd / defaultMetricValsColSqd
estBiasesPCSqd = (defaultBiasesApproxPC - defaultBiasesCol)**2
fracErrorPC = estBiasesPCSqd / defaultBiasesColSqd
#fracErrorPC = estBiasesPCSqd / defaultMetricValsColSqd


fracErrorMatrix = np.dstack((fracError, fracErrorPC)).squeeze()
df = pd.DataFrame(fracErrorMatrix,
                  index=metricsNames,
                  columns= ['fracError', 'fracErrorPC'])
errorsFig = px.line(df, x=df.index, y=df.columns,
              title = 'Fractional error in corrected solution.')
errorsFig.update_yaxes(title="Corr Err Sqd / Def Err Sqd")
errorsFig.update_xaxes(title="Metric and region")
errorsFig.update_layout(hovermode="x")


# Plot each column of normalized sensitivity matrix as a separate line.
# Each column tells us how all metrics vary with a single parameter.
df = pd.DataFrame(normlzdSensMatrix,
                  index=metricsNames,
                  columns=paramsNames)
normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized sensitivity matrix.<br>
                       Each column (line) shows how sensitive the metrics are to a change in a single parameter value.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""" )
normlzdSensMatrixColsFig.update_yaxes(title="Norml sens, (|param|/(obsmetric-defmetric)) * dmetric/dparam")
normlzdSensMatrixColsFig.update_xaxes(title="Metric and region")
normlzdSensMatrixColsFig.update_layout(hovermode="x")

# Plot each row of normalized sensitivity matrix as a separate line.
# Each row tells us how a single metric varies with all parameters.
df = pd.DataFrame(np.transpose(normlzdSensMatrix),
                  index=paramsNames,
                  columns=metricsNames)
normlzdSensMatrixRowsFig = px.line(df, x=df.index, y=df.columns,
              title = """Rows of normalized sensitivity matrix.<br>
                       Each row (line) tells us how a single metric varies with all parameters.<br>
                       (A positive value means that an increase in parameter value brings the default simulation closer to obs.)""" )
normlzdSensMatrixRowsFig.update_yaxes(title="Norml sens, (|param|/(obsmetric-defmetric)) * dmetric/dparam")
normlzdSensMatrixRowsFig.update_xaxes(title="Parameter")
normlzdSensMatrixRowsFig.update_layout(hovermode="x")

# Plot the parameter values recommended by SVD.
paramsMatrix = np.dstack((dparamsSoln,paramsSoln)).squeeze()
fracParamsMatrix = np.diagflat(np.reciprocal(np.abs(defaultParamValsOrigRow))) @ paramsMatrix
df = pd.DataFrame(fracParamsMatrix,
                  index=paramsNames,
                  columns= ['fracDparamsSoln', 'fracParamsSoln'])
fracParamsFig = px.line(df, x=df.index, y=df.columns,
              title = 'Fractional parameter values (and values - default) recommended by SVD.')
fracParamsFig.update_yaxes(title="(Parameter value) / abs(default value)")
fracParamsFig.update_xaxes(title="Parameter Name")
fracParamsFig.update_layout(hovermode="x")

# Plot the parameter values recommended by SVD.
# Multiply in the user-designated scale factors before plotting.
paramsFig = go.Figure()
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowVals[:,0]*paramsScales, name=r'$paramsSoln - \sigma$'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiVals[:,0]*paramsScales, fill='tonexty', name=r'$paramsSoln + \sigma$'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSoln[:,0]*paramsScales, name='paramsSoln'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales, name='paramsSolnPC'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales, name='defaultParamVals'))
paramsFig.update_yaxes(title="User-scaled parameter value")
paramsFig.update_xaxes(title="Parameter Name")
paramsFig.update_layout(hovermode="x")


sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),

        dcc.Graph( id='biasesFig', figure=biasesFig ),
        dcc.Graph( id='errorsFig', figure=errorsFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
        dcc.Graph( id='fracParamsFig', figure=fracParamsFig ),
        dcc.Graph( id='paramsFig', figure=paramsFig )

])


if __name__ == '__main__':
        sensMatrixDashboard.run_server(debug=True)