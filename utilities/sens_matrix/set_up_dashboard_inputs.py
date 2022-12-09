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
            analyzeSensMatrix, setupObsCol, setupDefaultMetricValsCol, \
            findOutliers, findParamsUsingElastic
    from test_analyzeSensMatrix import write_test_netcdf_files

    # L1 regularization coefficient, i.e., penalty on param perturbations in objFnc
    reglrCoef = 0.0

    # Metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsNamesAndWeights = [ \
                        ['SWCF_GLB', 8.00], \
                        ['SWCF_DYCOMS', 1.00], \
                        ['SWCF_HAWAII', 1.00], \
                        ['SWCF_VOCAL', 4.00], \
                        ['SWCF_VOCAL_near', 4.00], \
                        ['SWCF_LBA', 1.00], \
                        ['SWCF_WP', 1.00], \
#                        ['SWCF_EP', 1.00], \
#                        ['SWCF_NP', 1.00], \
#                        ['SWCF_SP', 1.00],  \
##                        ['SWCF_PA', 1.01], \
#                        ['SWCF_CAF', 1.00], \
                        ['SWCF_Nambian', 1.00], \
                        ['SWCF_Nambian_near', 1.00], \
                        ['LWCF_GLB', 8.00], \
###                        ['LWCF_DYCOMS', 1.01], \
###                        ['LWCF_HAWAII', 1.01], \
###                        ['LWCF_VOCAL', 1.01], \
#                        ['LWCF_LBA', 1.00], \
#                        ['LWCF_WP', 1.00], \
#                        ['LWCF_EP', 1.00], \
#                        ['LWCF_NP', 1.00], \
#                        ['LWCF_SP', 1.00], \
####                        ['LWCF_PA',  1.01], \
###                        ['LWCF_CAF', 1.01], \
                        ['PRECT_GLB', 8.00], \
#                        ['PRECT_LBA', 1.00], \
#                        ['PRECT_WP', 1.00] \
###                        ['PRECT_EP', 1.01], \
###                        ['PRECT_NP', 1.01], \
###                        ['PRECT_SP', 1.01], \
####                        ['PRECT_PA', 1.01], \
##                        ['PRECT_CAF', 1.00] \
                         ]

#                        ['PRECT_DYCOMS', 0.01], \
#                        ['PRECT_HAWAII', 0.01], \
#                        ['PRECT_VOCAL', 0.01], \

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
    folder_name = 'Regional_files/20221207_cam/'
    paramsNamesScalesAndFilenames = [ \
#                    ['clubb_c7', 1.0, \
#                     folder_name + 'zt_2_Regional.nc',  \
#                     folder_name + 'zt_3_Regional.nc'], \
#                    ['clubb_c11', 1.0, \
#                     folder_name + 'zt_c11m_Regional.nc',  \
#                     folder_name + 'zt_c11p_Regional.nc'], \
#                    ['clubb_gamma_coef', 1.0, \
#                     folder_name + 'zt_gamma_coefm_Regional.nc',  \
#                     folder_name + 'zt_gamma_coefp_Regional.nc'], \
#                    ['clubb_beta', 1.0, \
#                     folder_name + 'zt_betam_Regional.nc',  \
#                     folder_name + 'zt_betap_Regional.nc'], \
                    ['clubb_c8', 1.0, \
                     folder_name + 'zt_c8m_Regional.nc',  \
                     folder_name + 'zt_c8p_Regional.nc'], \
#                    ['clubb_c_uu_shr', 1.0, \
#                     folder_name + 'zt_c_uu_shrm_Regional.nc', \
#                     folder_name + 'zt_c_uu_shrp_Regional.nc'], \
#                    ['clubb_c_uu_buoy', 1.0, \
#                     folder_name + 'zt_c_uu_buoym_Regional.nc', \
#                     folder_name + 'zt_c_uu_buoyp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_n2m_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_n2p_Regional.nc'], \
#                    ['clubb_altitude_threshold', 0.001, \
#                     folder_name + 'zt_altitude_thresholdm_Regional.nc',
#                     folder_name + 'zt_altitude_thresholdp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_bkgnd', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_bkgndm_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_bkgndp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_sfc', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_sfcm_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_sfcp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_shear', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_shearm_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_shearp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3, \
#                     folder_name + 'zt_c_invrs_tau_wpxp_n2_threshm_Regional.nc', \
#                     folder_name + 'zt_c_invrs_tau_wpxp_n2_threshp_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2_wp2', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_n2_wp2m_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_n2_wp2p_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2_xp2', 1.0, \
#                         folder_name + 'zt_c_invrs_tau_n2_xp2m_Regional.nc',
#                         folder_name + 'zt_c_invrs_tau_n2_xp2p_Regional.nc'], \
#                    ['clubb_c_invrs_tau_n2_clear_wp3', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_n2_clear_wp3m_Regional.nc',
#                     folder_name + 'zt_c_invrs_tau_n2_clear_wp3p_Regional.nc'], \
#                    ['clubb_c_invrs_tau_wpxp_ri', 1.0, \
#                     folder_name + 'zt_c_invrs_tau_wpxp_rim_Regional.nc', \
#                     folder_name + 'zt_c_invrs_tau_wpxp_rip_Regional.nc'], \
                    ['micro_mg_autocon_fact', 1.0, \
                     folder_name + 'zt_autocon_factm_Regional.nc', \
                     folder_name + 'zt_autocon_factp_Regional.nc'], \
                    ['micro_mg_dcs', 1.0, \
                         folder_name + 'zt_dcsm_Regional.nc', \
                         folder_name + 'zt_dcsp_Regional.nc'], \
                        ]

    dfparamsNamesScalesAndFilenames =  \
        pd.DataFrame( paramsNamesScalesAndFilenames, \
                          columns = ['paramsNames', 'paramsScales',
                                     #'sensNcFilenames', 'sensNcFilenamesExt'] )
                                     'sensNcFilenamesExt', 'sensNcFilenames'] )
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
        folder_name + 'zt0_Regional.nc'
#        '20220903/anvil.bmg20220630.sens723_1_Regional.nc'

    # Metrics from simulation that use the SVD-recommended parameter values
    # Here, we use default simulation just as a placeholder.
    linSolnNcFilename = \
            folder_name + 'zt0_Regional.nc'
#            folder_name + 'zt_23_Regional.nc'

# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# NOTE: PRECT is in the unit of m/s
    obsMetricValsDict = { \
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
    'LWCF_Nambian': 12.3294, 'PRECT_Nambian':0.00000000177636 , 'SWCF_Nambian': -66.9495, 'TMQ_Nambian': 24.4823, \
    'LWCF_Nambian_near': 10.904, 'PRECT_Nambian_near':0.00000000238369 , 'SWCF_Nambian_near': -36.1216, 'TMQ_Nambian_near': 17.5188, \
        }

    # Set up a column vector of observed metrics
    obsMetricValsCol = setupObsCol(obsMetricValsDict, metricsNames)

    return (metricsNames, metricsWeights, \
        paramsNames, paramsScales, \
        transformedParamsNames, \
        sensNcFilenames, sensNcFilenamesExt, \
        defaultNcFilename, linSolnNcFilename, \
        obsMetricValsDict, obsMetricValsCol, \
        reglrCoef)


#if __name__ == '__main__':
#    main()
#        sensMatrixDashboard.run_server(debug=True)
