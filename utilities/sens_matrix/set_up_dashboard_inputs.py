# -*- coding: utf-8 -*-

# Run this app with `python3 sens_matrix_dashboard.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

"""
Set up input data to sens_matrix_dashboard.
This includes assigning variables for input netcdf filenames,
and setting regional metric weights and observed values of parameters.
"""


def setUpInputs():

#    import dash
#    import dash_core_components as dcc
#    import dash_html_components as html
#    import plotly.express as px
#    import plotly.graph_objects as go
    import pandas as pd

    import numpy as np
    import pdb
    import sklearn
#    from plotly.figure_factory import create_quiver
    from itertools import chain

    from analyze_sensitivity_matrix import \
            analyzeSensMatrix, \
            findOutliers, findParamsUsingElastic
    from test_analyzeSensMatrix import write_test_netcdf_files

    # L1 regularization coefficient, i.e., penalty on param perturbations in objFnc
    # Increase this value to 0.1 or 0.5 or so if you want to reduce 
    # the size of tuned parameter perturbations.
    reglrCoef = 0.0

    # Metrics are observed quantities that we want a tuned simulation to match.
    #    The first column is the metric name.
    #    The order of metricNames determines the order of rows in sensMatrix.
    # The second column is a vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fitting process.
    #   Use a large weight for global (GLB) metrics.
    # The third column is a vector of normalization values for metrics.  
    #   If a value in the 3rd column is set to -999, then the metric is simply normalized by the observed value.
    #   Otherwise, the value in the 3rd column is itself the normalization value for the metric.  
    metricsNamesWeightsAndNorms = [ \
#                        ['TMQ_RMSE', 1.00, 15.], \
#                        ['PSL_RMSE', 1.00, 1000.], \
#                        ['TS_RMSE', 1.00, 15.], \
#                        ['LHFLX_RMSE', 1.00, 15.], \
#                        ['SHFLX_RMSE', 1.00, 15.], \
#                        ['CLDLOW_RMSE', 1.00, 15.], \
#                        ['SWCF_RACC', 0.01, 0.2], \
#                        ['SWCF_RMSEP', 1.0, 15.], \
#                        ['SWCF_RMSE', 0.01, 15.], \
#                        ['SWCF_GLB', 1.00, -999], \
#                        ['SWCF_DYCOMS', 1.00, -999], \
#                        ['SWCF_HAWAII', 2.00, -999], \
#                        ['SWCF_VOCAL', 2.00, -999], \
#                        ['SWCF_VOCAL_near', 0.1, -999], \
#                        ['SWCF_LBA', 1.00, -999], \
#                        ['SWCF_WP', 1.00, -999], \
#                        ['SWCF_EP', 1.00, -999], \
#                        ['SWCF_NP', 1.00, -999], \
#                        ['SWCF_SP', 1.00, -999],  \
##                        ['SWCF_PA', 1.01, -999], \
#                        ['SWCF_CAF', 1.00, -999], \
#                        ['SWCF_Namibia', 1.00, -999], \
#                        ['SWCF_Namibia_near', 0.1, -999], \
    ['SWCF_1_1', 0.0056, -999], ['SWCF_1_2', 0.0056, -999],  ['SWCF_1_3', 0.0056, -999], ['SWCF_1_4', 0.0056, -999], \
    ['SWCF_1_5', 0.0056, -999], ['SWCF_1_6', 0.0056, -999],  ['SWCF_1_7', 0.0056, -999], ['SWCF_1_8', 0.0056, -999], \
    ['SWCF_1_9', 0.0056, -999], ['SWCF_1_10', 0.0056, -999], ['SWCF_1_11', 0.0056, -999], ['SWCF_1_12', 0.0056, -999], \
    ['SWCF_2_1', 0.015, -999], ['SWCF_2_2', 0.015, -999],  ['SWCF_2_3', 0.015, -999], ['SWCF_2_4', 0.015, -999], \
    ['SWCF_2_5', 0.015, -999], ['SWCF_2_6', 0.015, -999],  ['SWCF_2_7', 0.015, -999], ['SWCF_2_8', 0.015, -999], \
    ['SWCF_2_9', 0.015, -999], ['SWCF_2_10', 0.015, -999], ['SWCF_2_11', 0.015, -999], ['SWCF_2_12', 0.015, -999], \
    ['SWCF_3_1', 0.0208, -999], ['SWCF_3_2', 0.0208, -999],  ['SWCF_3_3', 0.0208, -999], ['SWCF_3_4', 0.0208, -999], \
    ['SWCF_3_5', 0.0208, -999], ['SWCF_3_6', 0.0208, -999],  ['SWCF_3_7', 0.0208, -999], ['SWCF_3_8', 0.0208, -999], \
    ['SWCF_3_9', 0.0208, -999], ['SWCF_3_10', 0.0208, -999], ['SWCF_3_11', 0.0208, -999], ['SWCF_3_12', 0.0208, -999], \
    ['SWCF_4_1', 0.0208, -999], ['SWCF_4_2', 0.0208, -999],  ['SWCF_4_3', 0.0208, -999], ['SWCF_4_4', 0.0208, -999], \
    ['SWCF_4_5', 0.0208, -999], ['SWCF_4_6', 0.0208, -999],  ['SWCF_4_7', 0.0208, -999], ['SWCF_4_8', 0.0208, -999], \
    ['SWCF_4_9', 0.0208, -999], ['SWCF_4_10', 0.0208, -999], ['SWCF_4_11', 0.0208, -999], ['SWCF_4_12', 0.0208, -999], \
    ['SWCF_5_1', 0.015, -999], ['SWCF_5_2', 0.015, -999],  ['SWCF_5_3', 0.015, -999], ['SWCF_5_4', 0.015, -999], \
    ['SWCF_5_5', 0.015, -999], ['SWCF_5_6', 0.015, -999],  ['SWCF_5_7', 0.015, -999], ['SWCF_5_8', 0.015, -999], \
    ['SWCF_5_9', 0.015, -999], ['SWCF_5_10', 0.015, -999], ['SWCF_5_11', 0.015, -999], ['SWCF_5_12', 0.015, -999], \
    ['SWCF_6_1', 0.0056, -999], ['SWCF_6_2', 0.0056, -999],  ['SWCF_6_3', 0.0056, -999], ['SWCF_6_4', 0.0056, -999], \
    ['SWCF_6_5', 0.0056, -999], ['SWCF_6_6', 0.0056, -999],  ['SWCF_6_7', 0.0056, -999], ['SWCF_6_8', 0.0056, -999], \
    ['SWCF_6_9', 0.0056, -999], ['SWCF_6_10', 0.0056, -999], ['SWCF_6_11', 0.0056, -999], ['SWCF_6_12', 0.0056, -999], \
#     ['SWCF_1_1', 0.41, -999], ['SWCF_1_2', 0.41, -999], ['SWCF_1_3', 0.41, -999], ['SWCF_1_4', 0.41, -999], \
#     ['SWCF_1_5', 0.41, -999], ['SWCF_1_6', 0.41, -999], ['SWCF_1_7', 0.41, -999], ['SWCF_1_8', 0.41, -999], \
#     ['SWCF_2_1', 1.00, -999], ['SWCF_2_2', 1.00, -999], ['SWCF_2_3', 1.00, -999], ['SWCF_2_4', 1.00, -999], \
#     ['SWCF_2_5', 1.00, -999], ['SWCF_2_6', 1.00, -999], ['SWCF_2_7', 1.00, -999], ['SWCF_2_8', 1.00, -999], \
#     ['SWCF_3_1', 1.00, -999], ['SWCF_3_2', 1.00, -999], ['SWCF_3_3', 1.00, -999], ['SWCF_3_4', 1.00, -999], \
#     ['SWCF_3_5', 1.00, -999], ['SWCF_3_6', 1.00, -999], ['SWCF_3_7', 1.00, -999], ['SWCF_3_8', 1.00, -999], \
#     ['SWCF_4_1', 0.41, -999], ['SWCF_4_2', 0.41, -999], ['SWCF_4_3', 0.41, -999], ['SWCF_4_4', 0.41, -999], \
#     ['SWCF_4_5', 0.41, -999], ['SWCF_4_6', 0.41, -999], ['SWCF_4_7', 0.41, -999], ['SWCF_4_8', 0.41, -999], \
#                        ['LWCF_GLB', 1.00, -999], \
###                        ['LWCF_DYCOMS', 1.01, -999], \
###                        ['LWCF_HAWAII', 1.01, -999], \
###                        ['LWCF_VOCAL', 1.01, -999], \
##                        ['LWCF_LBA', 1.00, -999], \
##                       ['LWCF_WP', 1.00, -999], \
###                        ['LWCF_EP', 1.01, -999], \
##                        ['LWCF_NP', 1.01, -999], \
##                        ['LWCF_SP', 1.01, -999], \
####                        ['LWCF_PA',  1.01, -999], \
###                        ['LWCF_CAF', 1.01, -999], \
#                        ['PRECT_GLB', 3.00, -999], \
#                        ['PRECT_RACC', 0.01, 1.0], \
#                        ['PRECT_RMSEP', 0.01, 1.0], \
#                        ['PRECT_RMSE', 0.01, 1.0], \
##                        ['PRECT_LBA', 1.00, -999], \
##                        ['PRECT_WP', 1.00, -999], \
###                        ['PRECT_EP', 1.01, -999], \
###                        ['PRECT_NP', 1.01, -999], \
###                        ['PRECT_SP', 1.01, -999], \
####                        ['PRECT_PA', 1.01, -999], \
##                        ['PRECT_CAF', 1.00, -999], \
#                        ['PSL_DYCOMS', 1.e0, 1e3], \
#                        ['PSL_HAWAII', 1.e0, 1e3], \
#                        ['PSL_VOCAL', 1.e0, 1e3], \
#                        ['PSL_VOCAL_near', 1.00, 1e3], \
#                        ['PSL_LBA', 1.e0, 1e3], \
#                        ['PSL_WP', 1.e0, 1e3], \
#                        ['PSL_EP', 1.e0, 1e3], \
#                        ['PSL_NP', 1.e0, 1e3], \
#                        ['PSL_SP', 1.e0, 1e3],  \
#                        ['PSL_PA', 1.00, 1e3], \
#                        ['PSL_CAF', 1.e0, 1e3], \
##                        ['PSL_Namibia', 1.00, 1e3], \
##                        ['PSL_Namibia_near', 1.00, 1e3], \
                         ]

#                        ['PRECT_DYCOMS', 0.01, -999], \
#                        ['PRECT_HAWAII', 0.01, -999], \
#                        ['PRECT_VOCAL', 0.01, -999], \

    # Split up the list above into metric names and the corresponding weights.
    dfMetricsNamesWeightsAndNorms =  \
        pd.DataFrame( metricsNamesWeightsAndNorms, columns = ['metricsNames', 'metricsWeights', 'metricsNorms'] )
    metricsNames = dfMetricsNamesWeightsAndNorms[['metricsNames']].to_numpy().astype(str)[:,0]
    metricsWeights = dfMetricsNamesWeightsAndNorms[['metricsWeights']].to_numpy().astype(float)
    metricsNorms = dfMetricsNamesWeightsAndNorms[['metricsNorms']].to_numpy().astype(float)


    # Parameters are tunable model parameters, e.g. clubb_C8.
    # The float listed below after the parameter name is a factor that is used below for scaling plots.
    #   It is not a weight and doesn't affect optimized values; it just makes the plots more readable.
    # Each parameter is associated with two sensitivity simulations; in one, the parameter is perturbed
    #    up and in the other, it is perturbed down.
    #    The output from each sensitivity simulation is expected to be stored in its own netcdf file.
    #    Each netcdf file contains metric values and parameter values for a single simulation.
    folder_name = 'Regional_files/20231204_30x30regs/'  # folder where regional netcdf files are stored.
    #folder_name = 'Regional_files/20231027/'  # folder where regional netcdf files are stored.
    #folder_name = 'Regional_files/20230910/'  # folder where regional netcdf files are stored.
    paramsNamesScalesAndFilenames = [ \
##                    ['clubb_c7', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_2.ne30pg2_r05_oECv3_Regional.nc',  \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_3.ne30pg2_r05_oECv3_Regional.nc'], \
##                    ['clubb_c11', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_4.ne30pg2_r05_oECv3_Regional.nc',  \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_5.ne30pg2_r05_oECv3_Regional.nc'], \
##                    ['clubb_gamma_coef', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_10.ne30pg2_r05_oECv3_Regional.nc',  \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_11.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c8', 1.0, \
                     folder_name + 'sens0707_14_Regional.nc',  \
                     folder_name + 'sens0707_15_Regional.nc'], \
                    ['clubb_c_k10', 1.0, \
                     folder_name + 'sens0707_12_Regional.nc', \
                     folder_name + 'sens0707_13_Regional.nc'], \
                    ['clubb_c_invrs_tau_n2', 1.0, \
                     folder_name + 'sens0707_10_Regional.nc',
                     folder_name + 'sens0707_11_Regional.nc'], \
                    ['clubb_altitude_threshold', 0.001, \
                     folder_name + 'sens0707_20_Regional.nc',
                     folder_name + 'sens0707_21_Regional.nc'], \
                    ['clubb_c_invrs_tau_sfc', 1.0, \
                     folder_name + 'sens0707_6_Regional.nc',
                     folder_name + 'sens0707_7_Regional.nc'], \
                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, \
                     folder_name + 'sens0707_8_Regional.nc', \
                     folder_name + 'sens0707_9_Regional.nc'], \
                    ['clubb_c_invrs_tau_n2_wp2', 1.0, \
                     folder_name + 'sens0707_4_Regional.nc',
                     folder_name + 'sens0707_5_Regional.nc'], \
##                    ['clubb_c_invrs_tau_wpxp_ri', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_16.ne30pg2_r05_oECv3_Regional.nc', \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_17.ne30pg2_r05_oECv3_Regional.nc'], \
##                    ['clubb_c_k1', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_25.ne30pg2_r05_oECv3_Regional.nc', \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_26.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_shear', 1.0, \
                     folder_name + 'sens0707_2_Regional.nc', \
                     folder_name + 'sens0707_3_Regional.nc'], \
##                    ['micro_vqit', 1.0, \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_16.ne30pg2_r05_oECv3_Regional.nc', \
##                     folder_name + 'chrysalis.bmg20220630.sens1107_17.ne30pg2_r05_oECv3_Regional.nc'], \
##                    ['clubb_c_invrs_tau_n2_xp2', 1.0, \
##                     folder_name + 'sens0707_18_Regional.nc', \
##                     folder_name + 'sens0707_19_Regional.nc'], \
##                     '20220903/anvil.bmg20220630.sens723_12.ne30pg2_r05_oECv3_Regional.nc',
##                     '20220903/anvil.bmg20220630.sens723_13.ne30pg2_r05_oECv3_Regional.nc'], \
                    ['clubb_c_invrs_tau_bkgnd', 1.0, \
                     folder_name + 'sens0707_16_Regional.nc',
                     folder_name + 'sens0707_17_Regional.nc'], \
##                    ['clubb_c_wp2_splat', 1.0, \
##                     folder_name + 'sens0707_26_Regional.nc',
##                     folder_name + 'sens0707_27_Regional.nc'], \
#                    ['cldfrc_dp2', 1.0, \
#                     folder_name + 'sens0707_38_Regional.nc',
#                     folder_name + 'sens0707_39_Regional.nc'], \
                        ]

    # Split up the above list into parameter names, scales, and filenames.
    dfparamsNamesScalesAndFilenames =  \
        pd.DataFrame( paramsNamesScalesAndFilenames, \
                          columns = ['paramsNames', 'paramsScales',
                                     'sensNcFilenames', 'sensNcFilenamesExt'] )
                                     #'sensNcFilenamesExt', 'sensNcFilenames'] )
    paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:,0]
    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:,0]
    sensNcFilenames = dfparamsNamesScalesAndFilenames[['sensNcFilenames']].to_numpy().astype(str)[:,0]
    sensNcFilenamesExt = dfparamsNamesScalesAndFilenames[['sensNcFilenamesExt']].to_numpy().astype(str)[:,0]

    # Below we designate the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    #transformedParamsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2', 'clubb_c_invrs_tau_n2_clear_wp3'])
    transformedParamsNames = np.array([''])

    prescribedParamsNamesScalesAndValues = \
                [ \
#                    ['clubb_c8', 1.0, 0.7, \
#                     folder_name + 'sens0707_14_Regional.nc',  \
#                     folder_name + 'sens0707_15_Regional.nc'], \
#                    ['clubb_c_k10', 1.0, 0.3, \
#                     folder_name + 'sens0707_12_Regional.nc', \
#                     folder_name + 'sens0707_13_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2', 1.0, 0.2, \
#                     folder_name + 'sens0707_10_Regional.nc',
#                     folder_name + 'sens0707_11_Regional.nc'], \
#                    ['clubb_c_invrs_tau_sfc', 1.0, 0.05, \
#                     folder_name + 'sens0707_6_Regional.nc',
#                     folder_name + 'sens0707_7_Regional.nc'], \
#                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, 0.00042, \
#                     folder_name + 'sens0707_8_Regional.nc', \
#                     folder_name + 'sens0707_9_Regional.nc'], \
#                    ['clubb_c_invrs_tau_shear', 1.0, 0.22, \
#                     folder_name + 'sens0707_2_Regional.nc', \
#                     folder_name + 'sens0707_3_Regional.nc'], \
#                    ['clubb_c_invrs_tau_bkgnd', 1.0, 1.1, \
#                     folder_name + 'sens0707_16_Regional.nc',
#                     folder_name + 'sens0707_17_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2_wp2', 1.0, 0.1, \
#                     folder_name + 'sens0707_4_Regional.nc',
#                     folder_name + 'sens0707_5_Regional.nc'], \
                ]
    # Split up the above list into parameter names, scales, and filenames.
    dfprescribedParamsNamesScalesAndValues =  \
        pd.DataFrame( prescribedParamsNamesScalesAndValues, \
                          columns = ['prescribedParamsNames', 
                                     'prescribedParamsScales',
                                     'prescribedParamVals',
                                     'prescribedSensNcFilenames', 'prescribedSensNcFilenamesExt'
                                    ] \
                    )
    prescribedParamsNames = dfprescribedParamsNamesScalesAndValues[['prescribedParamsNames']].to_numpy().astype(str)[:,0]
    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    prescribedParamsScales = dfprescribedParamsNamesScalesAndValues[['prescribedParamsScales']].to_numpy().astype(float)[:,0]
    prescribedParamVals = dfprescribedParamsNamesScalesAndValues[['prescribedParamVals']].to_numpy().astype(float)[:,0]
    prescribedParamValsRow = prescribedParamVals
    prescribedSensNcFilenames = dfprescribedParamsNamesScalesAndValues[['prescribedSensNcFilenames']].to_numpy().astype(str)[:,0]
    prescribedSensNcFilenamesExt = dfprescribedParamsNamesScalesAndValues[['prescribedSensNcFilenamesExt']].to_numpy().astype(str)[:,0]
    prescribedTransformedParamsNames = np.array([''])


    # Netcdf file containing metric and parameter values from the default simulation
    defaultNcFilename = \
        folder_name + 'sens0707_1_Regional.nc'
#        folder_name + 'chrysalis.bmg20220630.sens1107_1.ne30pg2_r05_oECv3_Regional.nc'
#        '20220903/anvil.bmg20220630.sens723_1.ne30pg2_r05_oECv3_Regional.nc'

    # Metrics from the global simulation that use the tuner-recommended parameter values
    linSolnNcFilename = \
           folder_name + 'sens0707_1_Regional.nc'
           #folder_name + 'sens0707_25_Regional.nc'
           #folder_name + 'sens0707_29_Regional.nc'
           # folder_name + 'chrysalis.bmg20220630.sens1107_30.ne30pg2_r05_oECv3_Regional.nc'
#            folder_name + 'chrysalis.bmg20220630.sens1107_23.ne30pg2_r05_oECv3_Regional.nc'

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    # These observed metrics will be matched as closely as possible by analyzeSensMatrix.
    # NOTE: PRECT is in the unit of m/s
    obsMetricValsDict = { \
    'SWCF_RACC': 0, \
    'SWCF_RMSEP': 0, \
    'SWCF_RMSE': 0, 'TMQ_RMSE': 0, 'PSL_RMSE': 0, 'TS_RMSE': 0, 'LHFLX_RMSE': 0, 'SHFLX_RMSE': 0, 'CLDLOW_RMSE': 0, \
    'LWCF_GLB': 28.008, 'PRECT_GLB': 0.000000031134259, 'SWCF_GLB': -45.81, 'TMQ_GLB': 24.423, \
    'LWCF_DYCOMS': 19.36681938, 'PRECT_DYCOMS':0.000000007141516, 'SWCF_DYCOMS': -63.49394226, 'TMQ_DYCOMS':20.33586884,\
    'LWCF_LBA': 43.83245087, 'PRECT_LBA':0.000000063727875, 'SWCF_LBA': -55.10041809, 'TMQ_LBA': 44.27890396,\
    'LWCF_HAWAII': 23.6855, 'PRECT_HAWAII':0.00000002087774, 'SWCF_HAWAII': -33.1536, 'TMQ_HAWAII': 32.4904,\
    'LWCF_WP': 54.5056, 'PRECT_WP':0.000000077433568, 'SWCF_WP': -62.3644, 'TMQ_WP':50.5412,\
    'LWCF_EP': 33.42149734, 'PRECT_EP': 0.000000055586694, 'SWCF_EP': -51.79394531, 'TMQ_EP':44.34251404,\
    'LWCF_NP': 26.23941231, 'PRECT_NP':0.000000028597503, 'SWCF_NP': -50.92364502, 'TMQ_NP':12.72111988,\
    'LWCF_SP': 31.96141052, 'PRECT_SP':0.000000034625369, 'SWCF_SP': -70.26461792, 'TMQ_SP':10.95032024,\
    'LWCF_PA': 47.32126999, 'PRECT_PA':0.000000075492694, 'SWCF_PA': -78.27433014, 'TMQ_PA':47.25967789,\
    'LWCF_CAF': 43.99757003784179687500, 'PRECT_CAF':0.000000042313699, 'SWCF_CAF': -52.50243378, 'TMQ_CAF':36.79592514,\
    'LWCF_VOCAL': 43.99757004, 'PRECT_VOCAL':0.000000001785546, 'SWCF_VOCAL': -77.26232147, 'TMQ_VOCAL':17.59922791, \
    'LWCF_VOCAL_near': 15.4783, 'PRECT_VOCAL_near':0.0000000037719, 'SWCF_VOCAL_near': -58.4732, 'TMQ_VOCAL_near': 14.9315, \
    'LWCF_Namibia': 12.3294, 'PRECT_Namibia':0.00000000177636 , 'SWCF_Namibia': -66.9495, 'TMQ_Namibia': 24.4823, \
    'LWCF_Namibia_near': 10.904, 'PRECT_Namibia_near':0.00000000238369 , 'SWCF_Namibia_near': -36.1216, 'TMQ_Namibia_near': 17.5188, \
    'PRECT_RACC': 0, \
    'PRECT_RMSEP': 0, \
    'PRECT_RMSE': 0, \
    'PSL_DYCOMS': 101868.515625, \
    'PSL_HAWAII': 101656.578125, \
    'PSL_VOCAL': 101668.703125, \
    'PSL_VOCAL_near': 101766.8203125, \
    'PSL_Namibia_near': 101741.7265625, \
    'PSL_Namibia far': 101550.6640625, \
    'PSL_LBA': 101052.40625, \
    'PSL_WP':  100909.4140625, \
    'PSL_EP':  101116.875, \
    'PSL_SP':  100021.4921875, \
    'PSL_NP':  101314.546875, \
    'PSL_PA':  100990.25, \
    'PSL_CAF': 100941.7890625, \
    'SWCF_1_1': -40.59888, 'SWCF_1_2': -38.27685, 'SWCF_1_3': -29.92496, 'SWCF_1_4': -28.98982, \
    'SWCF_1_5': -26.83143, 'SWCF_1_6': -25.31917, 'SWCF_1_7': -27.9917, 'SWCF_1_8': -24.87963, \
    'SWCF_1_9': -17.89297, 'SWCF_1_10': -18.54498, 'SWCF_1_11': -19.26345, 'SWCF_1_12': -37.08834, \
    'SWCF_2_1': -34.4901, 'SWCF_2_2': -30.40586, 'SWCF_2_3': -30.86093, 'SWCF_2_4': -39.89378, \
    'SWCF_2_5': -52.95258, 'SWCF_2_6': -71.01636, 'SWCF_2_7': -67.66965, 'SWCF_2_8': -58.65926, \
    'SWCF_2_9': -35.28358, 'SWCF_2_10': -49.01012, 'SWCF_2_11': -57.6655, 'SWCF_2_12': -47.73322, \
    'SWCF_3_1': -24.91453, 'SWCF_3_2': -18.33076, 'SWCF_3_3': -37.69804, 'SWCF_3_4': -64.86975, \
    'SWCF_3_5': -54.89531, 'SWCF_3_6': -42.75589, 'SWCF_3_7': -38.91958, 'SWCF_3_8': -53.20009, \
    'SWCF_3_9': -47.8519, 'SWCF_3_10': -49.38494, 'SWCF_3_11': -35.26117, 'SWCF_3_12': -31.96173, \
    'SWCF_4_1': -46.085, 'SWCF_4_2': -36.71035, 'SWCF_4_3': -43.04711, 'SWCF_4_4': -45.04954, \
    'SWCF_4_5': -40.85117, 'SWCF_4_6': -52.4172, 'SWCF_4_7': -44.30988, 'SWCF_4_8': -32.32601, \
    'SWCF_4_9': -36.84876, 'SWCF_4_10': -62.21399, 'SWCF_4_11': -46.25127, 'SWCF_4_12': -35.85495, \
    'SWCF_5_1': -64.23532, 'SWCF_5_2': -66.20067, 'SWCF_5_3': -66.24201, 'SWCF_5_4': -63.59531, \
    'SWCF_5_5': -55.99519, 'SWCF_5_6': -57.86568, 'SWCF_5_7': -60.74911, 'SWCF_5_8': -66.81461, \
    'SWCF_5_9': -63.63784, 'SWCF_5_10': -53.90995, 'SWCF_5_11': -60.69987, 'SWCF_5_12': -70.45469, \
    'SWCF_6_1': -24.44015, 'SWCF_6_2': -22.021, 'SWCF_6_3': -22.83143, 'SWCF_6_4': -19.53066, \
    'SWCF_6_5': -22.34129, 'SWCF_6_6': -30.43135, 'SWCF_6_7': -39.39627, 'SWCF_6_8': -37.68586, \
    'SWCF_6_9': -38.25595, 'SWCF_6_10': -31.7473, 'SWCF_6_11': -27.47777, 'SWCF_6_12': -27.86843
#    'SWCF_1_1': -42.18101, 'SWCF_1_2': -32.30865, 'SWCF_1_3': -32.55965, 'SWCF_1_4': -45.65011, \
#    'SWCF_1_5': -50.19252, 'SWCF_1_6': -32.94496, 'SWCF_1_7': -37.16001, 'SWCF_1_8': -47.75633, \
#    'SWCF_2_1': -24.65771, 'SWCF_2_2': -29.43244, 'SWCF_2_3': -58.16404, 'SWCF_2_4': -52.84988, \
#    'SWCF_2_5': -48.13301, 'SWCF_2_6': -48.1224, 'SWCF_2_7': -46.68948, 'SWCF_2_8': -35.85141, \
#    'SWCF_3_1': -46.92494, 'SWCF_3_2': -46.37671, 'SWCF_3_3': -46.52961, 'SWCF_3_4': -49.72062, \
#    'SWCF_3_5': -46.60182,  'SWCF_3_6': -42.66064, 'SWCF_3_7': -55.72744, 'SWCF_3_8': -45.17018, \
#    'SWCF_4_1': -52.65081, 'SWCF_4_2': -50.30189, 'SWCF_4_3': -46.82841, 'SWCF_4_4': -48.10809, \
#    'SWCF_4_5': -53.8271, 'SWCF_4_6': -54.2933, 'SWCF_4_7': -47.07664, 'SWCF_4_8': -54.32018 
    }

    # Comment out if not using 20x20reg files
    metricsNamesWeightsAndNorms = setUp_x_MetricsList(defaultNcFilename)
    # Split up the list above into metric names and the corresponding weights.
    dfMetricsNamesWeightsAndNorms =  \
        pd.DataFrame( metricsNamesWeightsAndNorms, columns = ['metricsNames', 'metricsWeights', 'metricsNorms'] )
    metricsNames = dfMetricsNamesWeightsAndNorms[['metricsNames']].to_numpy().astype(str)[:,0]
    metricsWeights = dfMetricsNamesWeightsAndNorms[['metricsWeights']].to_numpy().astype(float)
    metricsNorms = dfMetricsNamesWeightsAndNorms[['metricsNorms']].to_numpy().astype(float)
    obsMetricValsDict = setUp_x_ObsMetricValsDict(folder_name + "OBS.nc")

    return (metricsNames, metricsWeights, metricsNorms, \
            obsMetricValsDict, \
            paramsNames, paramsScales, \
            transformedParamsNames, \
            prescribedParamsNames, prescribedParamsScales, \
            prescribedTransformedParamsNames, \
            prescribedParamValsRow, \
            prescribedSensNcFilenames, prescribedSensNcFilenamesExt, \
            sensNcFilenames, sensNcFilenamesExt, \
            defaultNcFilename, linSolnNcFilename, \
            reglrCoef)


def setUpPreliminaries(metricsNames, metricsNorms, \
                       obsMetricValsDict, \
                       paramsNames, transformedParamsNames, \
                       prescribedParamsNames, prescribedParamValsRow, \
                       prescribedTransformedParamsNames, \
                       sensNcFilenames, \
                       defaultNcFilename \
                      ):

    import numpy as np
    import pdb
    import netCDF4

    # Set up a column vector of observed metrics
    obsMetricValsCol = setUpObsCol(obsMetricValsDict, metricsNames)

    # Set up a normalization vector for metrics, normMetricValsCol.
    # It equals the observed value when metricsNorms has the special value of -999, 
    #     but otherwise it is set manually in metricsNorms itself.
    normMetricValsCol = metricsNorms
    for idx in np.arange(len(metricsNorms)):
        if np.isclose(metricsNorms[idx],-999.0): 
            normMetricValsCol[idx] = obsMetricValsCol[idx]

    # Based on the default simulation,
    #    set up a row vector of parameter values.
    numParams = len(paramsNames)
    defaultParamValsRow, defaultParamValsOrigRow = \
            setupDefaultParamVectors(paramsNames, transformedParamsNames,
                                numParams,
                                defaultNcFilename)

    # Create row vector size numParams containing
    # parameter values from sensitivity simulations
    sensParamValsOrigRow = np.zeros((1, numParams))
    # This variable contains transformed parameter values,
    #    if transformedParamsNames is non-empty:
    sensParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        # Read netcdf file with changed parameter values from all sensitivity simulations.
        f_sensParams = netCDF4.Dataset(sensNcFilenames[idx], 'r')
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        sensParamValsOrigRow[0,idx] = f_sensParams.variables[paramName][0]
        # Transform [0,1] variable to extend over range [0,infinity]
        if paramName in transformedParamsNames:
            #sensParamValsRow[0,idx] = -np.log(1-sensParamValsRow[0,idx])
            sensParamValsRow[0,idx] = np.log(sensParamValsOrigRow[0,idx])
        else:
            sensParamValsRow[0,idx] = sensParamValsOrigRow[0,idx]
        f_sensParams.close()


    # Calculate the magnitude of the maximum value of parameters
    #    from the default run (and sensitivity runs as a backup), for later use
    #    in scaling the normalized sensitivity matrix.
    # Initially, set values to the default-simulation values
    magParamValsRow = np.abs(defaultParamValsRow)
    # Now replace any zero default values with the value from the sensitivity run
    for idx, elem in np.ndenumerate(defaultParamValsRow):
        if (np.abs(elem) <= np.finfo(elem.dtype).eps): # if default value is zero
            magParamValsRow[0,idx[1]] = np.abs(sensParamValsRow[0,idx[1]]) # set to sensitivity value
    if np.any( np.isclose(magParamValsRow, np.zeros((1,numParams))) ):
        print("\nsensParamValsRow =")
        print(sensParamValsRow)
        print("\nmagParamValsRow =")
        print(magParamValsRow)
        sys.exit("Error: A parameter value from both default and sensitivity simulation is zero.")

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    #print("defaultMetricValsCol=", defaultMetricValsCol)

    # Store biases in default simulation
    # defaultBiasesCol = + delta_b
    defaultBiasesCol = np.subtract(defaultMetricValsCol, obsMetricValsCol)

    #pdb.set_trace()

    # Based on the default simulation,
    #    set up a row vector of prescribed parameter values.
    numPrescribedParams = len(prescribedParamsNames)
    defaultPrescribedParamValsRow, defaultPrescribedParamValsOrigRow = \
            setupDefaultParamVectors(prescribedParamsNames, prescribedTransformedParamsNames,
                                     numPrescribedParams,
                                     defaultNcFilename)

    # Calculate the magnitude of the maximum value of parameters
    #    from the default run (and sensitivity runs as a backup), for later use
    #    in scaling the normalized sensitivity matrix.
    # Initially, set values to the default-simulation values
    magPrescribedParamValsRow = np.abs(defaultPrescribedParamValsRow)
    # Now replace any zero default values with the value from the sensitivity run
    for idx, elem in np.ndenumerate(defaultPrescribedParamValsRow):
        if (np.abs(elem) <= np.finfo(elem.dtype).eps): # if default value is zero
            magPrescribedParamValsRow[0,idx[1]] = np.abs(prescribedParamValsRow[0,idx[1]]) # set to prescribed value
    if np.any( np.isclose(magPrescribedParamValsRow, np.zeros((1,numPrescribedParams))) ):
        print("\nprescribedParamValsRow =")
        print(prescribedParamValsRow)
        print("\nmagPrescribedParamValsRow =")
        print(magPrescribedParamValsRow)
        sys.exit("Error: A prescribed parameter value is zero and so is the prescribed default value.")

    dnormlzdPrescribedParams = ( prescribedParamValsRow - defaultPrescribedParamValsRow ) \
                                / magPrescribedParamValsRow

    #print("prescribedParamValsRow=", prescribedParamValsRow)
    #print("defaultPrescribedParamValsRow=", defaultPrescribedParamValsRow)
    #print("magPrescribedParamValsRow=", magPrescribedParamValsRow)

    #pdb.set_trace()

    dnormlzdPrescribedParams = dnormlzdPrescribedParams.T


    return ( obsMetricValsCol, normMetricValsCol, \
             defaultBiasesCol, \
             defaultParamValsOrigRow, \
             magParamValsRow, \
             dnormlzdPrescribedParams, \
             magPrescribedParamValsRow
           )


def setUp_x_ObsMetricValsDict(obsPathAndFilename):
    """ 
    This is intended for the case in which 20x20deg fields are used.
    Input: Filename containing observed values of metrics.
    Output: Dictionary of observations.
    """

    import numpy as np
    import netCDF4
    import re

    # Read netcdf file with metrics and parameters from default simulation
    f_obs = netCDF4.Dataset(obsPathAndFilename, 'r')

    obsMetricValsDict = {}
    for varName in f_obs.variables:
        #print(varName)
        if re.search("^SWCF_[0-9]+_",varName) and not "MSWCF" in varName:
            varEntry = f_obs[varName]
            varVal = varEntry[:].data[:][0]
            obsMetricValsDict[varName] = varVal
            #print((varName, varVal))

    f_obs.close()

    #print(obsMetricValsDict)
    #print(len(obsMetricValsDict))

    return obsMetricValsDict


def setUp_x_MetricsList(defPathAndFilename):
    """
    This is intended for the case in which 20x20deg fields are used.
    Input: Filename of default run.
    Output: List of 20x20reg metric values.
    """

    import numpy as np
    import netCDF4
    import re

    # Read netcdf file with metrics and parameters from default simulation
    f_def = netCDF4.Dataset(defPathAndFilename, 'r')

    metricsNamesWeightsAndNorms = []
    for varName in f_def.variables:
        #print(varName)
        if re.match("^numb_[0-9]+_[0-9]+",varName):
            areaWeightEntry = f_def[varName]
            areaWeightVal = areaWeightEntry[:].data[:][0]
            SWCF_string = varName.replace("numb", "SWCF")
            metricsNamesWeightsAndNorms.append([SWCF_string,  areaWeightVal, -999])
            #print((SWCF_string, areaWeightVal))

    f_def.close()

    #print(obsMetricValsDict)
    #print(metricsNamesWeightsAndNorms)

    return metricsNamesWeightsAndNorms


def setUpObsCol(obsMetricValsDict, metricsNames):
    """
    Input: A python dictionary of observed metrics.
    Output: A column vector of observed metrics
    """

    import numpy as np
    import pdb 

    # Number of metrics
    numMetrics = len(metricsNames)

    # Set up column vector of numMetrics elements containing
    # "true" metric values from observations
    obsMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        obsMetricValsCol[idx] = obsMetricValsDict[metricName]

    return obsMetricValsCol


def setupDefaultParamVectors(paramsNames, transformedParamsNames,
                        numParams,
                        defaultNcFilename):
    """
    Input: Filename containing default-simulation metrics and parameters.
    Output: Row vector of default-sim parameter values.
    """

    import numpy as np
    import netCDF4

    # Read netcdf file with metrics and parameters from default simulation
    f_defaultMetricsParams = netCDF4.Dataset(defaultNcFilename, 'r')

    # Create row vector size numParams containing
    # parameter values from default simulation
    defaultParamValsOrigRow = np.zeros((1, numParams))
    defaultParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        defaultParamValsOrigRow[0,idx] = f_defaultMetricsParams.variables[paramName][0]
        # Transform [0,1] variable to extend over range [0,infinity]
        if paramName in transformedParamsNames:
            #defaultParamValsRow[0,idx] = -np.log(1-defaultParamValsOrigRow[0,idx])
            defaultParamValsRow[0,idx] = np.log(defaultParamValsOrigRow[0,idx])
        else:
            defaultParamValsRow[0,idx] = defaultParamValsOrigRow[0,idx]

    f_defaultMetricsParams.close()

    return (defaultParamValsRow, defaultParamValsOrigRow)

def setupDefaultMetricValsCol(metricsNames, defaultNcFilename):
    """
    Input: Filename containing default-simulation metrics.
    Output: Column vector of default-sim metrics.
    """

    import numpy as np
    import netCDF4

    # Number of metrics
    numMetrics = len(metricsNames)

    # Read netcdf file with metrics and parameters from default simulation
    f_defaultMetricsParams = netCDF4.Dataset(defaultNcFilename, 'r')

    # Set up column vector of numMetrics elements containing
    # metric values from default simulation
    defaultMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        defaultMetricValsCol[idx] = f_defaultMetricsParams.variables[metricName][0]

    f_defaultMetricsParams.close()

    return defaultMetricValsCol



#if __name__ == '__main__':
#    main()
#        sensMatrixDashboard.run_server(debug=True)
