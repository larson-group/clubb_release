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

from analyze_sensitivity_matrix import analyzeSensMatrix, setupObsCol, setupDefaultMetricValsCol
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
                        ['SWCF_GLB', 2.01], \
                        ['SWCF_DYCOMS', 1.01], \
                        ['SWCF_HAWAII', 1.01], \
                        ['SWCF_VOCAL', 1.01], \
                        ['SWCF_LBA', 0.01], \
                        ['SWCF_WP', 4.01], \
                        ['SWCF_EP', 4.01], \
                        ['SWCF_NP', 0.01], \
                        ['SWCF_SP', 0.01], \
                        ['SWCF_PA', 0.01], \
                        ['SWCF_CAF', 0.01], \
                        ['LWCF_GLB', 0.01], \
#                        ['LWCF_DYCOMS', 0.01], \
#                        ['LWCF_HAWAII', 0.01], \
#                        ['LWCF_VOCAL', 0.01], \
                        ['LWCF_LBA', 0.01], \
                        ['LWCF_WP', 0.01], \
                        ['LWCF_EP', 0.01], \
                        ['LWCF_NP', 0.01], \
                        ['LWCF_SP', 0.01], \
                        ['LWCF_PA', 0.01], \
                        ['LWCF_CAF', 0.01], \
                        ['PRECT_GLB', 2.01], \
#                        ['PRECT_DYCOMS', 0.01], \
#                        ['PRECT_HAWAII', 0.01], \
#                        ['PRECT_VOCAL', 0.01], \
                        ['PRECT_LBA', 0.01], \
                        ['PRECT_WP', 0.01], \
                        ['PRECT_EP', 0.01], \
                        ['PRECT_NP', 0.01], \
                        ['PRECT_SP', 0.01], \
                        ['PRECT_PA', 0.01], \
                        ['PRECT_CAF', 0.01] \
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
#                  ['clubb_c8', 1.0, '20210515/anvil.newp3.vqitp4_c81.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210515/anvil.newp3.vqitp4_thres2p8.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_n2', 1.0, '20210515/anvil.newp3.vqitp4_n2p5.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_k10', 1.0, '20210515/anvil.newp3.vqitp4_ck10p3.ne30_ne30_Regional.nc'], \
#                  ['vqit', 1.0, '20210515/anvil.newp3.vqitp35.ne30_ne30_Regional.nc'], \
#                  ['clubb_c_invrs_tau_wpxp_ri', 1.0, '20210515/anvil.newp3.vqitp4_ri5.ne30_ne30_Regional.nc'], \
#                  ['clubb_gamma_coef', 1.0, '20210515/anvil.newp3.vqitp4_gap2.ne30_ne30_Regional.nc'], \
#                  ['max_total_ni', 1.0, '20210515/anvil.newp3.vqitp4_ni700.ne30_ne30_Regional.nc'], \
#                  ['cldfrc_rhminl', 1.0, '20210515/anvil.newp3.vqitp4_rhp9.ne30_ne30_Regional.nc'] \
#                    ['clubb_c8', 1.0, '20210705/anvil.0703.c8p3.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210705/anvil.0703.n2thresp3.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['clubb_c_invrs_tau_n2', 1.0, '20210705/anvil.0703.n2p2.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['clubb_c_k10', 1.0, '20210705/anvil.0703.ck10p7.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['vqit', 1.0, '20210705/anvil.0703.vqitp5.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c8', 1.0, '20210706/anvil.0703.tuner0706_c8p44.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210706/anvil.0703.tuner0706_n2thres3p2.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_n2', 1.0, '20210706/anvil.0703.tuner0706_n2p3.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_k10', 1.0, '20210706/anvil.0703.tuner0706_ck10p32.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['vqit', 1.0, '20210706/anvil.0703.tuner0706_vqitp5.ne30pg2_r05_oECv3_Regional.nc'], \
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
        '20210706/anvil.0703.tuner0706.ne30pg2_r05_oECv3_Regional.nc'
       # '20210705/anvil.0703.newdefault_n2thresp25.ne30pg2_r05_oECv3_Regional.nc'

# Metrics from simulation that use the SVD-recommended parameter values
# Here, we use default simulation just as a placeholder.
linSolnNcFilename = \
        '20210706/anvil.0703.tuner0706.ne30pg2_r05_oECv3_Regional.nc'
       # '20210705/anvil.0703.newdefault_n2thresp25.ne30pg2_r05_oECv3_Regional.nc'
       #'20210705/anvil.0703.tuner0706.ne30pg2_r05_oECv3_Regional.nc'


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
defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
normlzdWeightedDefaultBiasesApproxPC, \
defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
defaultParamValsOrigRow, dparamsSoln, \
paramsSoln, paramsLowVals, paramsHiVals, \
paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
         analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

# Set up a column vector of observed metrics
obsMetricValsCol = setupObsCol(obsMetricValsDict, metricsNames)

# Set up a column vector of metric values from the default simulation
linSolnMetricValsCol = setupDefaultMetricValsCol(metricsNames, linSolnNcFilename)

# Set up a column vector of metric values from the default simulation
defaultMetricValsCol = setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

# Store biases in default simulation
linSolnBiasesCol = np.subtract(linSolnMetricValsCol, defaultMetricValsCol)

#weightedBias = metricsWeights * ( 1.0 - defaultBiasesApprox / defaultBiasesCol )
#weightedBiasMag = np.dot( np.transpose(weightedBias), weightedBias )
#weightedBiasMagRatio = weightedBiasMag / np.dot( np.transpose(metricsWeights), metricsWeights )
## weightedBiasPC = metricsWeights * ( obs - forward ) / ( obs - def ) ?
#weightedBiasPC = metricsWeights * ( 1.0 - defaultBiasesApproxPC / defaultBiasesCol )
#weightedBiasPCMag = np.dot( np.transpose(weightedBiasPC), weightedBiasPC )
#weightedBiasPCMagRatio = weightedBiasPCMag / np.dot( np.transpose(metricsWeights), metricsWeights )

fracBiasPC = ( defaultBiasesApproxPC - defaultBiasesCol ) * np.reciprocal(np.abs(defaultMetricValsCol))
fracDefaultBiasesCol = defaultBiasesCol * np.reciprocal(np.abs(defaultMetricValsCol))
fracBiasPCMagRatio = np.linalg.norm(fracBiasPC)**2 / np.linalg.norm(fracDefaultBiasesCol)**2

weightedBiasPC = normlzdWeightedDefaultBiasesApproxPC - metricsWeights
weightedBiasPCMagRatio = np.linalg.norm(weightedBiasPC)**2 / np.linalg.norm(metricsWeights)**2

weightedBiasLin = metricsWeights * ( 1.0 - linSolnBiasesCol / defaultBiasesCol )
weightedBiasLinMag = np.dot( np.transpose(weightedBiasLin), weightedBiasLin )
weightedBiasLinMagRatio = weightedBiasLinMag / np.dot( np.transpose(metricsWeights), metricsWeights )

#FracBiasPCMag = np.dot( np.transpose(defaultBiasesApproxPC), defaultBiasesApproxPC ) / \
#                np.dot( np.transpose(defaultBiasesCol), defaultBiasesCol )
weightedFracBiasPCMag = np.dot( np.transpose(metricsWeights*defaultBiasesApproxPC), metricsWeights*defaultBiasesApproxPC ) / \
                np.dot( np.transpose(metricsWeights*defaultBiasesCol), metricsWeights*defaultBiasesCol )

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Plot the biases of the default simulation and the SVD approximation of that
biasesMatrix = np.dstack((defaultBiasesCol,#defaultBiasesApprox,
                          defaultBiasesApproxPC,
                          defaultBiasesApproxLowValsPC,
                          defaultBiasesApproxHiValsPC,
                          linSolnBiasesCol)).squeeze()
fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(defaultMetricValsCol))) @ biasesMatrix
df = pd.DataFrame(fracBiasesMatrix,
                  index=metricsNames,
                  columns= ['fracDefBias',
                            #'fracDefBiasApprox',
                            'fracDefBiasesApproxPC',
                            'fracDefBiasesApproxLowValsPC',
                            'fracDefBiasesApproxHiValsPC',
                            'fracLinSolnBiasesCol'])
biasesFig = px.line(df, x=df.index, y=df.columns,
              title = """Fractional biases of default simulation and approximations thereof.<br>
                    Plotted quantities have the structure -(def-obs), -(def-fwd), -(def-lin)""")
biasesFig.update_yaxes(title="-(Def-Sim) / abs(default metric value)")
biasesFig.update_xaxes(title="Metric and region")
biasesFig.update_layout(hovermode="x")
#biasesFig.data[3].name = "fracDefBiasesApprox, " + "{:.2f}".format(weightedBiasMagRatio[0][0])
biasesFig.data[1].name = "fracDefBiasesApproxPC, " \
                        + "{:.2f}".format(weightedBiasPCMagRatio) \
                        + ", {:.2f}".format(fracBiasPCMagRatio)
biasesFig.data[4].name = "fracLinSolnBiasesCol, " + "{:.2f}".format(weightedBiasLinMagRatio[0][0])
#pdb.set_trace()


## Compute errors in the final parameter solution as a frac of the original errors
## estBiasesOrig = default metric value + correction - observed metric value
#defaultBiasesColSqd = defaultBiasesCol**2
##defaultMetricValsColSqd = defaultMetricValsCol**2
#estBiasesSqd = (defaultBiasesApprox - defaultBiasesCol)**2
#fracError = estBiasesSqd / defaultBiasesColSqd
##fracError = estBiasesSqd / defaultMetricValsColSqd
#estBiasesPCSqd = (defaultBiasesApproxPC - defaultBiasesCol)**2
#fracErrorPC = estBiasesPCSqd / defaultBiasesColSqd
##fracErrorPC = estBiasesPCSqd / defaultMetricValsColSqd


#fracErrorMatrix = np.dstack((fracError, fracErrorPC)).squeeze()
#df = pd.DataFrame(fracErrorMatrix,
#                  index=metricsNames,
#                   columns= ['fracError', 'fracErrorPC'])
#errorsFig = px.line(df, x=df.index, y=df.columns,
#              title = 'Fractional error in corrected solution.')
#errorsFig.update_yaxes(title="Corr Err Sqd / Def Err Sqd")
#errorsFig.update_xaxes(title="Metric and region")
#errorsFig.update_layout(hovermode="x")

# Plot the maximum with respect to each parameter of the normalized sensitivity of each metric.
df = pd.DataFrame(np.max(np.abs(normlzdSensMatrix), axis=1), # max of absolute val of each row
                  index=metricsNames,
                  columns=['Max abs normlzd sensitivity'])
maxSensMetricsFig = px.line(df, x=df.index, y=df.columns,
              title = """Maximum normalized sensitivity of each metric with respect to parameters.<br>
                       (Low sensitivity means that the metric is unbudgeable by these parameters.)""" )
maxSensMetricsFig.update_yaxes(title="Max |sens row|")
maxSensMetricsFig.update_xaxes(title="Metric and region")
maxSensMetricsFig.update_layout(hovermode="x")

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

## Plot the fractional parameter values recommended by SVD.
#paramsMatrix = np.dstack((dparamsSoln,paramsSoln)).squeeze()
#fracParamsMatrix = np.diagflat(np.reciprocal(np.abs(defaultParamValsOrigRow))) @ paramsMatrix
#df = pd.DataFrame(fracParamsMatrix,
#                  index=paramsNames,
#                  columns= ['fracDparamsSoln', 'fracParamsSoln'])
#fracParamsFig = px.line(df, x=df.index, y=df.columns,
#              title = 'Fractional parameter values (and values - default) recommended by SVD.')
#fracParamsFig.update_yaxes(title="(Parameter value) / abs(default value)")
#fracParamsFig.update_xaxes(title="Parameter Name")
#fracParamsFig.update_layout(hovermode="x")

# Plot the parameter values recommended by SVD.
# Multiply in the user-designated scale factors before plotting.
paramsFig = go.Figure()
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPC[:,0]*paramsScales,
                               name=r'$paramsSolnPC - \sigma$'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPC[:,0]*paramsScales, fill='tonexty',
                               name=r'$paramsSolnPC + \sigma$'))
#paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales, name='paramsSoln'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales, name='paramsSolnPC'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales, name='defaultParamVals'))
paramsFig.update_yaxes(title="User-scaled parameter value")
paramsFig.update_xaxes(title="Parameter Name")
paramsFig.update_layout(hovermode="x")


sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),

        dcc.Graph( id='biasesFig', figure=biasesFig ),
#        dcc.Graph( id='errorsFig', figure=errorsFig ),
        dcc.Graph( id='maxSensMetricsFig', figure=maxSensMetricsFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
#        dcc.Graph( id='fracParamsFig', figure=fracParamsFig ),
        dcc.Graph( id='paramsFig', figure=paramsFig )

])


if __name__ == '__main__':
        sensMatrixDashboard.run_server(debug=True)
