# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.


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
#                  ['clubb_c_invrs_tau_wpxp_ri', 1.0, '20210515/anvil.newp3.vqitp4_ri5.ne30_ne30_Regional.nc'], \
#                  ['clubb_gamma_coef', 1.0, '20210515/anvil.newp3.vqitp4_gap2.ne30_ne30_Regional.nc'], \
#                  ['max_total_ni', 1.0, '20210515/anvil.newp3.vqitp4_ni700.ne30_ne30_Regional.nc'], \
#                  ['cldfrc_rhminl', 1.0, '20210515/anvil.newp3.vqitp4_rhp9.ne30_ne30_Regional.nc'] \
#                   ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210705/anvil.0703.n2thresp3.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['clubb_c_invrs_tau_n2', 1.0, '20210705/anvil.0703.n2p2.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['clubb_c_k10', 1.0, '20210705/anvil.0703.ck10p7.ne30pg2_r05_oECv3_Regional.nc'], \
#                   ['vqit', 1.0, '20210705/anvil.0703.vqitp5.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c8', 1.0, '20210810/anvil.0703.tuner0706_2.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, '20210810/anvil.0703.tuner0706_5.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_n2', 1.0, '20210810/anvil.0703.tuner0706_4.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_k10', 1.0, '20210810/anvil.0703.tuner0706_3.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['micro_vqit', 1.0, '20210810/anvil.0703.tuner0706_6.ne30pg2_r05_oECv3_Regional.nc'], \
                        ]

dfparamsNamesScalesAndFilenames =  \
        pd.DataFrame( paramsNamesScalesAndFilenames, columns = ['paramsNames', 'paramsScales', 'sensNcFilenames'] )
paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:,0]
# Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
# The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:,0]
sensNcFilenames = dfparamsNamesScalesAndFilenames[['sensNcFilenames']].to_numpy().astype(str)[:,0]

# This the subset of paramsNames that vary from [0,1] (e.g., C5)
#    and hence will be transformed to [0,infinity] in order to make
#    the relationship between parameters and metrics more linear:
#transformedParamsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2', 'clubb_c_invrs_tau_n2_clear_wp3'])
transformedParamsNames = np.array([''])

# Netcdf file containing metric and parameter values from the default simulation
defaultNcFilename = \
        '20210810/anvil.0703.tuner0706_1.ne30pg2_r05_oECv3_Regional.nc'
       # '20210706/anvil.0703.tuner0706.ne30pg2_r05_oECv3_Regional.nc'
       # '20210705/anvil.0703.newdefault_n2thresp25.ne30pg2_r05_oECv3_Regional.nc'

# Metrics from simulation that use the SVD-recommended parameter values
# Here, we use default simulation just as a placeholder.
linSolnNcFilename = \
        '20210706/anvil.0703.tuner0706.ne30pg2_r05_oECv3_Regional.nc'
       # '20210705/anvil.0703.newdefault_n2thresp25.ne30pg2_r05_oECv3_Regional.nc'


# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# Global mean, NOTE PRECT is in the unit of m/s
#obsMetricValsDict = { \
#   'LWCF_GLB': 28.008, 'PRECT_GLB': 0.000000031134259, 'SWCF_GLB': -45.81, 'TMQ_GLB': 24.423, \
#   'LWCF_DYCOMS': 17.39229, 'PRECT_DYCOMS': 5.540433451401774e-09, 'SWCF_DYCOMS': -63.05767, 'TMQ_DYCOMS': 19.03481,\
#   'LWCF_LBA': 54.9896, 'PRECT_LBA': 8.75291e-08, 'SWCF_LBA': -67.0833, 'TMQ_LBA': 48.1414,\
#   'LWCF_HAWAII': 21.111, 'PRECT_HAWAII': 1.47959e-08, 'SWCF_HAWAII': -31.6793, 'TMQ_HAWAII': 30.7126,\
#   'LWCF_WP': 56.3396, 'PRECT_WP': 7.65602e-08, 'SWCF_WP': -64.7297, 'TMQ_WP':50.7996,\
#   'LWCF_EP': 32.724, 'PRECT_EP': 5.93702e-08, 'SWCF_EP': -54.0756, 'TMQ_EP':45.7974,\
#   'LWCF_NP': 30.3083, 'PRECT_NP': 3.657e-08, 'SWCF_NP': -73.3256, 'TMQ_NP':12.9728,\
#   'LWCF_SP': 30.6713, 'PRECT_SP': 3.11686e-08, 'SWCF_SP': -65.3633, 'TMQ_SP':12.5768,\
#   'LWCF_PA': 28.5308, 'PRECT_PA': 3.22142e-08, 'SWCF_PA': -73.4048, 'TMQ_PA':50.0178,\
#   'LWCF_CAF': 57.0782, 'PRECT_CAF': 5.71253e-08, 'SWCF_CAF': -60.8234, 'TMQ_CAF':41.1117,\
#   'LWCF_VOCAL': 16.2196, 'PRECT_VOCAL': 1.78555e-09, 'SWCF_VOCAL': -77.2623, 'TMQ_VOCAL':17.5992  }

obsMetricValsDict = { \
    'LWCF_GLB': 28.008, 'PRECT_GLB': 0.000000031134259, 'SWCF_GLB': -45.81, 'TMQ_GLB': 24.423, \
    'LWCF_DYCOMS': 17.39229, 'PRECT_DYCOMS':5.540433451401774e-09, 'SWCF_DYCOMS': -63.05767, 'TMQ_DYCOMS':19.03481,\
    'LWCF_LBA': 54.7332, 'PRECT_LBA':7.86887e-08, 'SWCF_LBA': -62.0982, 'TMQ_LBA': 51.1303,\
    'LWCF_HAWAII': 21.111, 'PRECT_HAWAII':1.47959e-08, 'SWCF_HAWAII': -31.6793, 'TMQ_HAWAII': 30.7126,\
    'LWCF_WP': 54.7332, 'PRECT_WP':7.86887e-08, 'SWCF_WP': -62.0982, 'TMQ_WP':51.1303,\
    'LWCF_EP': 32.724, 'PRECT_EP': 5.93702e-08, 'SWCF_EP': -54.0756, 'TMQ_EP':45.7974,\
    'LWCF_NP': 26.2394, 'PRECT_NP':2.85975e-08, 'SWCF_NP': -50.9236, 'TMQ_NP':12.7285,\
    'LWCF_SP': 31.9614, 'PRECT_SP':3.46254e-08, 'SWCF_SP': -70.2646, 'TMQ_SP':10.947,\
    'LWCF_PA': 28.5308, 'PRECT_PA':3.22142e-08, 'SWCF_PA': -73.4048, 'TMQ_PA':50.0178,\
    'LWCF_CAF': 57.0782, 'PRECT_CAF':5.71253e-08, 'SWCF_CAF': -60.8234, 'TMQ_CAF':41.1117,\
    'LWCF_VOCAL': 16.2196, 'PRECT_VOCAL':1.78555e-09, 'SWCF_VOCAL': -77.2623, 'TMQ_VOCAL':17.5992  }

# Calculate changes in parameter values needed to match metrics.
defaultMetricValsCol, defaultBiasesCol, \
defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
normlzdWeightedDefaultBiasesApprox, normlzdWeightedDefaultBiasesApproxPC, \
defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
vhNormlzd, uNormlzd, sNormlzd, \
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

# Calculate the fraction of the default-sim bias that is removed by tuning.
# This is unweighted and hence is not necessarily less than one.
# Bias = J*delta_p = ( fwd - obs ) = numerator
Bias = ( defaultBiasesApprox + defaultBiasesCol )
# DefaultBiasesCol = -delta_b = - ( default - obs ) = denominator
BiasMagRatio = np.linalg.norm(Bias)**2 / np.linalg.norm(defaultBiasesCol)**2

# Calculate the fraction of the default-sim bias that is removed by tuning.
# This is unweighted and hence is not necessarily less than one.
# fracBiasPC = J*delta_p = ( fwd - obs ) = numerator
BiasPC = ( defaultBiasesApproxPC + defaultBiasesCol )
# fracDefaultBiasesCol = -delta_b = - ( default - obs ) = denominator
#fracDefaultBiasesCol = defaultBiasesCol
BiasPCMagRatio = np.linalg.norm(BiasPC)**2 / np.linalg.norm(defaultBiasesCol)**2

# Calculate the fraction of bias removed, but normalized and weighted,
# like the equations that the SVD actually solves, so that according to theory,
# the value should be < 1.
# But I'm not sure if it will be < 1 if the parameters are transformed to log space.
normlzdMDeltaB = metricsWeights * defaultBiasesCol / np.abs(obsMetricValsCol) # right-hand side
weightedBiasNumer = normlzdWeightedDefaultBiasesApprox + normlzdMDeltaB
weightedBiasDenom = normlzdMDeltaB
weightedBiasMagRatio = np.linalg.norm(weightedBiasNumer)**2 / np.linalg.norm(weightedBiasDenom)**2

# Calculate the fraction of bias removed by PC soln, but normalized and weighted,
# like the equations that the SVD actually solves, so that according to theory,
# the value should be < 1.
# But I'm not sure if it will be < 1 if the parameters are transformed to log space.
weightedBiasPCNumer = normlzdWeightedDefaultBiasesApproxPC + normlzdMDeltaB
weightedBiasPCDenom = normlzdMDeltaB
weightedBiasPCMagRatio = np.linalg.norm(weightedBiasPCNumer)**2 / np.linalg.norm(weightedBiasPCDenom)**2

# Fraction of bias that is removed by the tuned, non-linear forward solution
# weightedBiasLin = metricsWeights * ( fwd - obs ) / (def - obs )
#weightedBiasLin = metricsWeights * ( 1.0 - linSolnBiasesCol / defaultBiasesCol )
#weightedBiasLinMag = np.dot( np.transpose(weightedBiasLin), weightedBiasLin )
#weightedBiasLinMagRatio = weightedBiasLinMag / np.dot( np.transpose(metricsWeights), metricsWeights )
# weightedBiasLinMag = numerator
#weightedBiasLinMag = np.dot( np.transpose(weightedBiasLin), weightedBiasLin )
#weightedBiasLinMagRatio = weightedBiasLinMag / \
#                          np.dot( np.transpose(metricsWeights*defaultBiasesCol),
#                                  metricsWeights*defaultBiasesCol )

# weightedBiasLin = metricsWeights * ( lin - obs ) = numerator
#weightedBiasLin = metricsWeights * ( linSolnBiasesCol - defaultBiasesCol )
#weightedBiasLinMagRatio = np.linalg.norm(weightedBiasLin)**2 / np.linalg.norm(weightedDefaultBiasesCol)**2


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Plot the biases of the default simulation and the SVD approximation of that
biasesMatrix = np.dstack((-defaultBiasesCol,
                          defaultBiasesApprox,
                          defaultBiasesApproxPC #,
                          #defaultBiasesApproxLowValsPC,
                          #defaultBiasesApproxHiValsPC,
                          #linSolnBiasesCol
                         )).squeeze()
fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(obsMetricValsCol))) @ biasesMatrix
df = pd.DataFrame(fracBiasesMatrix,
                  index=metricsNames,
                  columns= ['fracDefBias',
                            'fracDefBiasesApprox',
                            'fracDefBiasesApproxPC',
                            #'fracDefBiasesApproxLowValsPC',
                            #'fracDefBiasesApproxHiValsPC',
                            #'fracLinSolnBiasesCol'
                           ])
biasesFig = px.line(df, x=df.index, y=df.columns,
              title = """Fractional biases of default simulation and approximations thereof.<br>
                    Plotted quantities have the structure -(def-obs), -(def-fwd), -(def-lin)""")
biasesFig.update_yaxes(title="-(Def-Sim) / abs(obs metric value)")
biasesFig.update_xaxes(title="Metric and region")
biasesFig.layout.legend.title = "Default or which approximation"
biasesFig.update_layout(hovermode="x")
biasesFig.data[1].name = "fracDefBiasesApprox, " \
                        + "{:.2f}".format(weightedBiasMagRatio) \
                        + ", {:.2f}".format(BiasMagRatio)
biasesFig.data[2].name = "fracDefBiasesApproxPC, " \
                        + "{:.2f}".format(weightedBiasPCMagRatio) \
                        + ", {:.2f}".format(BiasPCMagRatio)
#biasesFig.data[5].name = "fracLinSolnBiasesCol, " + "{:.2f}".format(weightedBiasLinMagRatio)
#pdb.set_trace()


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
normlzdSensMatrixColsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
normlzdSensMatrixColsFig.update_xaxes(title="Metric and region")
normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
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
normlzdSensMatrixRowsFig.update_yaxes(title="Norml sens, (|param|/|obsmetric|) * dmetric/dparam")
normlzdSensMatrixRowsFig.update_xaxes(title="Parameter")
normlzdSensMatrixRowsFig.layout.legend.title = "Metric"
normlzdSensMatrixRowsFig.update_layout(hovermode="x")

# Plot each column of right-singular vector matrix, V.
#pdb.set_trace()
rightSingVectorNums = (np.arange(paramsNames.shape[0])+1).astype(str)
df = pd.DataFrame(np.transpose(vhNormlzd),
                  index=paramsNames,
                  columns=rightSingVectorNums)
vhNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted right-singular vector matrix, V.<br>
                       Each column (line) is a vector of parameter values associated with a singular value.<br>""" )
vhNormlzdColsFig.update_yaxes(title="Right-singular vector")
vhNormlzdColsFig.update_xaxes(title="Parameter")
vhNormlzdColsFig.layout.legend.title = "Singular vector"
vhNormlzdColsFig.update_layout(hovermode="x")
for idx, val in np.ndenumerate(sNormlzd):
    vhNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

# Plot each column of left-singular vector matrix, U.
#pdb.set_trace()
#leftSingVectorNums = (np.arange(metricsNames.shape[0])+1).astype(str)
df = pd.DataFrame(uNormlzd,
                  index=metricsNames,
                  columns=rightSingVectorNums)
uNormlzdColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted left-singular vector matrix, U.<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
uNormlzdColsFig.update_yaxes(title="Left-singular vector")
uNormlzdColsFig.update_xaxes(title="Parameter")
uNormlzdColsFig.layout.legend.title = "Singular vector"
uNormlzdColsFig.update_layout(hovermode="x")
for idx, val in np.ndenumerate(sNormlzd):
    uNormlzdColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

# Plot each column of left-singular vector matrix, U.
#pdb.set_trace()
#leftSingVectorNums = (np.arange(metricsNames.shape[0])+1).astype(str)
df = pd.DataFrame(uNormlzd*defaultBiasesCol/np.abs(obsMetricValsCol),
                  index=metricsNames,
                  columns=rightSingVectorNums)
uNormlzdBiasColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, unweighted left-singular vector matrix, U, dotted with delta_b (defaultBiasesCol).<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
uNormlzdBiasColsFig.update_yaxes(title="uNormlzd dot delta_b")
uNormlzdBiasColsFig.update_xaxes(title="Parameter")
uNormlzdBiasColsFig.layout.legend.title = "Singular vector"
uNormlzdBiasColsFig.update_layout(hovermode="x")
for idx, val in np.ndenumerate(sNormlzd):
    uNormlzdBiasColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)
#pdb.set_trace()

# Plot the parameter values recommended by SVD.
# Multiply in the user-designated scale factors before plotting.
paramsFig = go.Figure()
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsLowValsPC[:,0]*paramsScales,
                               name=r'$paramsSolnPC - \sigma$'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsHiValsPC[:,0]*paramsScales, fill='tonexty',
                               name=r'$paramsSolnPC + \sigma$'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales, name='paramsSoln'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales, name='paramsSolnPC'))
paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales, name='defaultParamVals'))
paramsFig.update_yaxes(title="User-scaled parameter value")
paramsFig.update_xaxes(title="Parameter Name")
paramsFig.update_layout(hovermode="x")

sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),

        dcc.Graph( id='biasesFig', figure=biasesFig ),
        dcc.Graph( id='maxSensMetricsFig', figure=maxSensMetricsFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
        dcc.Graph( id='vhNormlzdColsFig', figure=vhNormlzdColsFig ),
        dcc.Graph( id='uNormlzdColsFig', figure=uNormlzdColsFig ),
        dcc.Graph( id='uNormlzdBiasColsFig', figure=uNormlzdBiasColsFig ),
        dcc.Graph( id='paramsFig', figure=paramsFig )

])


if __name__ == '__main__':
        sensMatrixDashboard.run_server(debug=True)
