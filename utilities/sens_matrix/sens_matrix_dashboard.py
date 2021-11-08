# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.

def main():

    import dash
    import dash_core_components as dcc
    import dash_html_components as html
    import plotly.express as px
    import plotly.graph_objects as go
    import pandas as pd

    import numpy as np
    import pdb
    import sklearn

    from analyze_sensitivity_matrix import \
            analyzeSensMatrix, setupObsCol, setupDefaultMetricValsCol, \
            findOutliers, findParamsUsingElastic
    from test_analyzeSensMatrix import write_test_netcdf_files

    # Metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsNamesAndWeights = [ \
                        ['SWCF_GLB', 1.01], \
                        ['SWCF_DYCOMS', 0.01], \
                        ['SWCF_HAWAII', 0.01], \
                        ['SWCF_VOCAL', 1.01], \
                        ['SWCF_LBA', 1.01], \
                        ['SWCF_WP', 1.01], \
                        ['SWCF_EP', 1.01], \
                        ['SWCF_NP', 1.01], \
                        ['SWCF_SP', 1.01], \
                        ['SWCF_PA', 0.01], \
                        ['SWCF_CAF', 1.01], \
                        ['LWCF_GLB', 1.01], \
                        ['LWCF_DYCOMS', 0.01], \
                        ['LWCF_HAWAII', 0.01], \
                        ['LWCF_VOCAL', 0.01], \
                        ['LWCF_LBA', 0.01], \
                        ['LWCF_WP', 1.01], \
                        ['LWCF_EP', 1.01], \
                        ['LWCF_NP', 1.01], \
                        ['LWCF_SP', 1.01], \
                        ['LWCF_PA', 0.01], \
                        ['LWCF_CAF', 0.01], \
                        ['PRECT_GLB', 1.01], \
#                        ['PRECT_DYCOMS', 0.01], \
#                        ['PRECT_HAWAII', 0.01], \
#                        ['PRECT_VOCAL', 0.01], \
                        ['PRECT_LBA', 0.01], \
                        ['PRECT_WP', 1.01], \
                        ['PRECT_EP', 1.01], \
                        ['PRECT_NP', 1.01], \
                        ['PRECT_SP', 1.01], \
#                        ['PRECT_PA', 1.01], \
#                        ['PRECT_CAF', 1.01] \
                         ]

    dfMetricsNamesAndWeights =  \
        pd.DataFrame( metricsNamesAndWeights, columns = ['metricsNames', 'metricsWeights'] )
    metricsNames = dfMetricsNamesAndWeights[['metricsNames']].to_numpy().astype(str)[:,0]
    metricsWeights = dfMetricsNamesAndWeights[['metricsWeights']].to_numpy()


    # Parameters are tunable model parameters, e.g. clubb_C8.
    # The float listed below is a factor that is used below for scaling plots.
    # Each parameter is associated with two sensitivity simulations in which that parameter is perturbed
    #    either up or down.
    #    The output from each sensitivity simulation is expected to be stored in its own netcdf file.
    #    Each netcdf file contains metric values and parameter values for a single simulation.
    paramsNamesScalesAndFilenames = [ \
#                    ['clubb_c7', 1.0, \
#                     '20211010/anvil.0703.lmm_2.ne30pg2_r05_oECv3_Regional.nc',  \
#                     '20211010/anvil.0703.lmm_3.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c11', 1.0, \
                     '20211010/anvil.0703.lmm_5.ne30pg2_r05_oECv3_Regional.nc',  \
                     '20211010/anvil.0703.lmm_4.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_gamma_coef', 1.0, \
                     '20211010/anvil.0703.lmm_6.ne30pg2_r05_oECv3_Regional.nc',  \
                     '20211010/anvil.0703.lmm_7.ne30pg2_r05_oECv3_Regional.nc'], \
#                    ['clubb_c8', 1.0, \
#                     '20211010/anvil.0703.lmm_9.ne30pg2_r05_oECv3_Regional.nc',  \
#                     '20211010/anvil.0703.lmm_8.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_k10', 1.0, \
                     '20211010/anvil.0703.lmm_10.ne30pg2_r05_oECv3_Regional.nc', \
                     '20211010/anvil.0703.lmm_11.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_n2', 1.0, \
                     '20211010/anvil.0703.lmm_12.ne30pg2_r05_oECv3_Regional.nc',
                     '20211010/anvil.0703.lmm_13.ne30pg2_r05_oECv3_Regional.nc'], \
#                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, \
#                     '20211010/anvil.0703.lmm_14.ne30pg2_r05_oECv3_Regional.nc', \
#                     '20211010/anvil.0703.lmm_15.ne30pg2_r05_oECv3_Regional.nc'], \
#                    ['micro_vqit', 1.0, \
#                     '20211010/anvil.0703.lmm_16.ne30pg2_r05_oECv3_Regional.nc', \
#                     '20211010/anvil.0703.lmm_17.ne30pg2_r05_oECv3_Regional.nc'], \
                        ]

    dfparamsNamesScalesAndFilenames =  \
        pd.DataFrame( paramsNamesScalesAndFilenames, \
                          columns = ['paramsNames', 'paramsScales',
                                     'sensNcFilenames', 'sensNcFilenamesExt'] )
    paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:,0]
    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:,0]
    sensNcFilenames = dfparamsNamesScalesAndFilenames[['sensNcFilenames']].to_numpy().astype(str)[:,0]
    sensNcFilenamesExt = dfparamsNamesScalesAndFilenames[['sensNcFilenamesExt']].to_numpy().astype(str)[:,0]

    # This the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    #transformedParamsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2', 'clubb_c_invrs_tau_n2_clear_wp3'])
    transformedParamsNames = np.array([''])

    # Netcdf file containing metric and parameter values from the default simulation
    defaultNcFilename = \
        '20211010/anvil.0703.lmm_1.ne30pg2_r05_oECv3_Regional.nc'

    # Metrics from simulation that use the SVD-recommended parameter values
    # Here, we use default simulation just as a placeholder.
    linSolnNcFilename = \
            '20211010/anvil.0703.lmm_20.ne30pg2_r05_oECv3_Regional.nc'


# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# NOTE: PRECT is in the unit of m/s
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


    # Estimate non-linearity of the global model to perturbations in parameter values.
    # To do so, calculate radius of curvature of the three points from the default simulation
    #   and the two sensitivity simulations.
    #calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames,
    #                      metricsWeights,
    #                      sensNcFilenames, sensNcFilenamesExt, defaultNcFilename)


    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, \
    defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
    normlzdWeightedDefaultBiasesApprox, normlzdWeightedDefaultBiasesApproxPC, \
    defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, \
    normlzdWeightedSensMatrix, biasNormlzdSensMatrix, svdInvrsNormlzdWeighted, \
    vhNormlzd, uNormlzd, sNormlzd, \
    vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted, \
    magParamValsRow, \
    defaultParamValsOrigRow, dparamsSoln, dnormlzdParamsSoln, \
    dparamsSolnPC, dnormlzdParamsSolnPC, \
    paramsSoln, paramsLowVals, paramsHiVals, \
    paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
          analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                            metricsWeights,
                            sensNcFilenames, defaultNcFilename,
                            obsMetricValsDict)


    # Set up a column vector of observed metrics
    obsMetricValsCol = setupObsCol(obsMetricValsDict, metricsNames)

    # Create scatterplot to look at outliers
    #createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, obsMetricValsCol, metricsNames, paramsNames)

    # Find outliers by use of the ransac algorithm
    outlier_mask, defaultBiasesApproxRansac, normlzdWeightedDefaultBiasesApproxRansac, \
    dnormlzdParamsSolnRansac, paramsSolnRansac = \
        findOutliers(normlzdSensMatrix, normlzdWeightedSensMatrix, \
                     defaultBiasesCol, obsMetricValsCol, magParamValsRow, defaultParamValsOrigRow)
    print( "ransac_outliers = ", metricsNames[outlier_mask] )
    print( "ransac_inliers = ", metricsNames[~outlier_mask] )
    #pdb.set_trace()

    # Find best-fit params by use of the Elastic Net algorithm
    defaultBiasesApproxElastic, normlzdWeightedDefaultBiasesApproxElastic, \
    dnormlzdParamsSolnElastic, paramsSolnElastic = \
        findParamsUsingElastic(normlzdSensMatrix, normlzdWeightedSensMatrix, \
                     defaultBiasesCol, obsMetricValsCol, metricsWeights, magParamValsRow, defaultParamValsOrigRow)
    #pdb.set_trace()

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Set up a column vector of metric values from the global simulation based on optimized
    #     parameter values.
    linSolnMetricValsCol = setupDefaultMetricValsCol(metricsNames, linSolnNcFilename)

    # Store biases in default simulation, ( global_model - default )
    linSolnBiasesCol = np.subtract(linSolnMetricValsCol, defaultMetricValsCol)

    # Calculate the fraction of the default-sim bias that remains after tuning.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    Bias = ( defaultBiasesApprox + defaultBiasesCol )
    # defaultBiasesCol = -delta_b = - ( default - obs ) = -denominator
    BiasMagRatio = np.linalg.norm(Bias/np.abs(obsMetricValsCol))**2 / \
                   np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxPC = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasPC = ( defaultBiasesApproxPC + defaultBiasesCol )
    # defaultBiasesCol = -delta_b = - ( default - obs ) = -denominator
    BiasPCMagRatio = np.linalg.norm(BiasPC/np.abs(obsMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxRansac = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasRansac = ( defaultBiasesApproxRansac + defaultBiasesCol )
    # defaultBiasesCol = -delta_b = - ( default - obs ) = -denominator
    BiasRansacMagRatio = np.linalg.norm(BiasRansac/np.abs(obsMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of the default-sim bias that remains after tuning,
    #    but using a truncated PC observation.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApproxElastic = J*delta_p = ( fwd - def )
    # numerator = ( fwd - def ) + ( def - obs ) = ( fwd - obs )
    BiasElastic = ( defaultBiasesApproxElastic + defaultBiasesCol )
    # defaultBiasesCol = -delta_b = - ( default - obs ) = -denominator
    BiasElasticMagRatio = np.linalg.norm(BiasElastic/np.abs(obsMetricValsCol))**2 / \
                     np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the global-model bias relative to the default-sim bias.
    # This is unweighted and hence is not necessarily less than one.
    # defaultBiasesApprox = J*delta_p = ( fwd - def )
    # numerator = ( linSoln - def ) + ( def - obs ) = ( linSoln - obs )
    linSolnBias = ( linSolnBiasesCol + defaultBiasesCol )
    # defaultBiasesCol = -delta_b = - ( default - obs ) = -denominator
    linSolnBiasMagRatio = np.linalg.norm(linSolnBias/np.abs(obsMetricValsCol))**2 / \
                          np.linalg.norm(defaultBiasesCol/np.abs(obsMetricValsCol))**2

    # Calculate the fraction of bias removed by the non-PC soln, but normalized and weighted,
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

    # Calculate the fraction of bias removed by Ransac soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    weightedBiasRansacNumer = normlzdWeightedDefaultBiasesApproxRansac + normlzdMDeltaB
    weightedBiasRansacDenom = normlzdMDeltaB
    weightedBiasRansacMagRatio = np.linalg.norm(weightedBiasRansacNumer)**2 / np.linalg.norm(weightedBiasRansacDenom)**2

    # Calculate the fraction of bias removed by Elastic soln, but normalized and weighted,
    # like the equations that the SVD actually solves, so that according to theory,
    # the value should be < 1.
    # But I'm not sure if it will be < 1 if the parameters are transformed to log space.
    weightedBiasElasticNumer = normlzdWeightedDefaultBiasesApproxElastic + normlzdMDeltaB
    weightedBiasElasticDenom = normlzdMDeltaB
    weightedBiasElasticMagRatio = np.linalg.norm(weightedBiasElasticNumer)**2 / np.linalg.norm(weightedBiasElasticDenom)**2

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
    weightedBiasLinSoln = metricsWeights * ( linSolnBiasesCol + defaultBiasesCol ) / np.abs(obsMetricValsCol)
    weightedBiasLinSolnMagRatio = np.linalg.norm(weightedBiasLinSoln)**2 / np.linalg.norm(normlzdMDeltaB)**2


    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    sensMatrixDashboard = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    # Plot the biases of the default simulation and the SVD approximation of that
    biasesMatrix = np.dstack((-defaultBiasesCol,
                          defaultBiasesApprox,
                          defaultBiasesApproxPC,
                          defaultBiasesApproxRansac,
                          defaultBiasesApproxElastic,
                          linSolnBiasesCol
                         )).squeeze()
    fracBiasesMatrix = np.diagflat(np.reciprocal(np.abs(obsMetricValsCol))) @ biasesMatrix
    df = pd.DataFrame(fracBiasesMatrix,
                  index=metricsNames,
                  columns= ['fracDefBias',
                            'fracDefBiasesApprox',
                            'fracDefBiasesApproxPC',
                            'fracDefBiasesApproxRansac',
                            'fracDefBiasesApproxElastic',
                            'fracLinSolnBiasesCol'
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
    biasesFig.data[3].name = "fracDefBiasesApproxRansac, " \
                        + "{:.2f}".format(weightedBiasRansacMagRatio) \
                         + ", {:.2f}".format(BiasRansacMagRatio)
    biasesFig.data[4].name = "fracDefBiasesApproxElastic, " \
                        + "{:.2f}".format(weightedBiasElasticMagRatio) \
                         + ", {:.2f}".format(BiasElasticMagRatio)
    biasesFig.data[5].name = "fracLinSolnBiasesCol, " \
                        + "{:.2f}".format(weightedBiasLinSolnMagRatio) \
                         + ", {:.2f}".format(linSolnBiasMagRatio)
    #pdb.set_trace()


    # Plot a scatterplot of default-simulation bias and SVD approximation of that bias.
    # Each column tells us how all metrics vary with a single parameter.
    biasSensMatrix = np.concatenate((-defaultBiasesCol/np.abs(obsMetricValsCol),
                                     defaultBiasesApproxPC/np.abs(obsMetricValsCol)), axis=1)
                                     #defaultBiasesApproxPC/np.abs(obsMetricValsCol), normlzdSensMatrix), axis=1)
    biasAndParamsNames = ["bias", "bias_approx_pc"]
    #biasAndParamsNames = np.append(["bias", "bias_approx_pc"], paramsNames)
    df = pd.DataFrame(biasSensMatrix,
                  index=metricsNames,
                  columns=biasAndParamsNames)
    biasSensMatrixScatter = px.scatter(df, x="bias_approx_pc", y="bias", text=metricsNames,
    #biasSensMatrixScatter = px.scatter(df, x=np.append(["bias_approx_pc"], paramsNames), y="bias",
              title = """Columns of normalized sensitivity matrix.<br>
                       vs. bias vector.<br>
                       """ )
    biasSensMatrixOneOneLine = px.line(df, x="bias", y="bias")
    biasSensMatrixScatterFig = go.Figure(data=biasSensMatrixScatter.data
                                              + biasSensMatrixOneOneLine.data)
    biasSensMatrixScatterFig.update_yaxes(title="bias")
    biasSensMatrixScatterFig.update_xaxes(title="normlzd sens matrix elements")
    biasSensMatrixScatterFig.update_traces(textposition='top center')
    #normlzdSensMatrixColsFig.layout.legend.title = "Parameter"
    #pdb.set_trace()

    # Plot each row of normalized sensitivity matrix as a separate line.
    # Each row tells us how a single metric varies with all parameters.
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
    maxSensMetricsFig.update_traces(mode='lines+markers')

    # Plot each column of normalized sensitivity matrix as a separate line.
    # Each column tells us how all metrics vary with a single parameter.
    df = pd.DataFrame(normlzdSensMatrix,
                  index=metricsNames,
                  columns=paramsNames)
    normlzdSensMatrixColsFig = px.line(df, x=df.index, y=df.columns,
              title =  """Columns of normalized, unweighted sensitivity matrix.<br>
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
              title = """Rows of normalized, unweighted sensitivity matrix.<br>
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

    # Plot each column of normalized, unweighted left-singular vector matrix, U.
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

    # Plot each column of normalized, weighted left-singular vector matrix, U.
    df = pd.DataFrame(uNormlzdWeighted,
                  index=metricsNames,
                  columns=rightSingVectorNums)
    uNormlzdWeightedColsFig = px.line(df, x=df.index, y=df.columns,
              title = """Columns of normalized, weighted left-singular vector matrix, U.<br>
                       Each column (line) is a vector of metric values associated with a singular value.<br>""" )
    uNormlzdWeightedColsFig.update_yaxes(title="Left-singular vector")
    uNormlzdWeightedColsFig.update_xaxes(title="Parameter")
    uNormlzdWeightedColsFig.layout.legend.title = "Singular vector"
    uNormlzdWeightedColsFig.update_layout(hovermode="x")
    for idx, val in np.ndenumerate(sNormlzdWeighted):
        uNormlzdWeightedColsFig.data[idx[0]].name = "{}".format(idx[0]+1) + ", " + "{:.2e}".format(val)

    # Plot each column of left-singular vector matrix, U, multiplied by biases.
    # Plot each column of left-singular vector matrix, U, multiplied by biases.
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
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSoln[:,0]*paramsScales,
                                   name='paramsSoln, |dp|=' + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSoln)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnPC[:,0]*paramsScales,
                                   name='paramsSolnPC, |dpPC|=' + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnPC)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnRansac[:,0]*paramsScales,
                                   name='paramsSolnRansac, |dpRansac|=' + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnRansac)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=paramsSolnElastic[:,0]*paramsScales,
                                   name='paramsSolnElastic, |dpElastic|='
                                   + '{:.2e}'.format(np.linalg.norm(dnormlzdParamsSolnElastic)) ))
    paramsFig.add_trace(go.Scatter(x=paramsNames, y=defaultParamValsOrigRow[0,:]*paramsScales, name='defaultParamVals'))
    paramsFig.update_yaxes(title="User-scaled parameter value")
    paramsFig.update_xaxes(title="Parameter Name")
    paramsFig.update_layout(hovermode="x")

    sensMatrixDashboard.layout = html.Div(children=[
        html.H1(children='Sensitivity matrix diagnostics'),

        html.Div(children=''' '''),

        dcc.Graph( id='biasesFig', figure=biasesFig ),
        dcc.Graph( id='biasesSensScatterFig', figure=biasSensMatrixScatterFig ),
        dcc.Graph( id='maxSensMetricsFig', figure=maxSensMetricsFig ),
        dcc.Graph( id='normlzdSensMatrixColsFig', figure=normlzdSensMatrixColsFig ),
        dcc.Graph( id='normlzdSensMatrixRowsFig', figure=normlzdSensMatrixRowsFig ),
        dcc.Graph( id='vhNormlzdColsFig', figure=vhNormlzdColsFig ),
        dcc.Graph( id='uNormlzdWeightedColsFig', figure=uNormlzdWeightedColsFig ),
        dcc.Graph( id='uNormlzdColsFig', figure=uNormlzdColsFig ),
        dcc.Graph( id='uNormlzdBiasColsFig', figure=uNormlzdBiasColsFig ),
         dcc.Graph( id='paramsFig', figure=paramsFig )

    ])

    sensMatrixDashboard.run_server(debug=True)

    return

def calcNormlzdRadiusCurv(metricsNames, paramsNames, transformedParamsNames,
                          metricsWeights,
                          sensNcFilenames, sensNcFilenamesExt, defaultNcFilename):

    """
    Calculate radius of curvature of output from 2 sensitivity simulations plus the default
    simulation.
    """
    import numpy as np
    import sys
    import netCDF4
    import matplotlib.pyplot as plt
    import pdb

    from analyze_sensitivity_matrix import setupDefaultMetricValsCol, setupDefaultParamVectors, \
                                           setupSensArrays

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Number of tunable parameters
    numParams = len(paramsNames)

    # Number of metrics
    numMetrics = len(metricsNames)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
        setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
                                 numMetrics, numParams,
                                 defaultNcFilename)

    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenames)

    # Set up sensitivity-simulation matrices from the extended sensitivity simulation
    sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
        setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        sensNcFilenamesExt)

    normlzd_radius_of_curv = np.full_like(sensMetricValsMatrix, 0.0)

    # Calculate differences in parameter values between default, sensitivity,
    #    and extended sensitivity runs.
    delta_params_def_sens = sensParamValsRow - defaultParamValsRow
    delta_params_def_sensExt = sensParamValsRowExt - defaultParamValsRow
    delta_params_sens_sensExt = sensParamValsRowExt - sensParamValsRow

    # Calculate numMetrics x numParams matrix of metric values.
    delta_metrics_def_sens = sensMetricValsMatrix - defaultMetricValsMatrix
    delta_metrics_def_sensExt = sensMetricValsMatrixExt - defaultMetricValsMatrix
    delta_metrics_sens_sensExt = sensMetricValsMatrixExt - sensMetricValsMatrix
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):
            # Distance between points in simulations = sqrt(dparam**2 + dmetric**2)
            length_def_sens = np.linalg.norm([delta_params_def_sens[0][col],
                                              delta_metrics_def_sens[row][col]])
            length_def_sensExt = np.linalg.norm([delta_params_def_sensExt[0][col],
                                                 delta_metrics_def_sensExt[row][col]])
            length_sens_sensExt = np.linalg.norm([delta_params_sens_sensExt[0][col],
                                                  delta_metrics_sens_sensExt[row][col]])
            semi_perim = 0.5 * ( length_def_sens + length_def_sensExt + length_sens_sensExt )
            # area of triangle formed by points.  Use Heron's formula.
            area = np.sqrt( semi_perim *
                           (semi_perim-length_def_sens) *
                           (semi_perim-length_def_sensExt) *
                           (semi_perim-length_sens_sensExt)
                          )
            if (area == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, area == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Greatest distance between parameter values in the 3 simulations:
            max_params_width = \
            np.max(np.abs([delta_params_def_sens[0][col],
                        delta_params_def_sensExt[0][col],
                        delta_params_sens_sensExt[0][col]]))
            if (max_params_width == 0.0):
                print( '\nIn calcNormlzdRadiusCurv, max_params_width == 0.0 for param ', paramsNames[col],
                        'and metric ', metricsNames[row] )

            # Calculate Menger curvature from triangle area and distance between points:
            normlzd_radius_of_curv[row][col] = 0.25 * length_def_sens*length_def_sensExt*length_sens_sensExt \
                                                / area / max_params_width



    fig, axs = plt.subplots(numMetrics, numParams)
    for col in np.arange(numParams):
        for row in np.arange(numMetrics):

            paramVals = [defaultParamValsRow[0][col], sensParamValsRow[0][col], sensParamValsRowExt[0][col]]
            metricVals = [defaultMetricValsMatrix[row][col], sensMetricValsMatrix[row][col],
                  sensMetricValsMatrixExt[row][col]]

            axs[row, col].scatter( paramVals, metricVals )
            axs[row, col].set_xlabel(paramsNames[col])
            axs[row, col].set_ylabel(metricsNames[row])
            fig.show()

    #pdb.set_trace()

    return

def createPcaBiplot(normlzdSensMatrix, defaultBiasesCol, obsMetricValsCol, metricsNames, paramsNames):

    import numpy as np
    from pca import pca
    import pdb

    # reduce the data towards 2 PCs
    model = pca(n_components=2, detect_outliers='ht2')

    # Augmented array with LHS and RHS
    augMatrix = np.concatenate((normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) ), axis=1)

    paramsList = list(paramsNames)
    paramsList.append('dbias')
    augParamsNames = np.asarray(paramsList)

    # Fit transform
    results = model.fit_transform(augMatrix, row_labels=metricsNames, col_labels=augParamsNames)
    #results = model.fit_transform(normlzdSensMatrix, row_labels=metricsNames, col_labels=paramsNames)

    PC_test = model.transform(augMatrix)
    #PC_test = model.transform(normlzdSensMatrix)
    outliers, outliers_params = model.compute_outliers(PC=PC_test)
    print("PCA outliers = ", outliers)

    # Plot explained variance
    #fig, ax = model.plot()

    # Scatter first 2 PCs
    #fig, ax = model.scatter()

    # Make biplot with the number of features
    fig, ax = model.biplot(n_feat=paramsNames.size+1)
    #fig, ax = model.biplot(n_feat=paramsNames.size)

    #pdb.set_trace()

    return

if __name__ == '__main__':
    main()
#        sensMatrixDashboard.run_server(debug=True)
