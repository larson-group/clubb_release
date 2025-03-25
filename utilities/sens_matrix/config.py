# -*- coding: utf-8 -*-

# Run this app with `python3 quadtune_driver.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

"""
In this file, users may specify input data to quadtune_driver.
This includes assigning filenames for input netcdf files,
regional metric weights, and observed values of parameters.
"""

import numpy as np
import pandas as pd


def setUpConfig(beVerbose):
    from set_up_inputs import (
        setUp_x_MetricsList,
        setupDefaultMetricValsCol, setUp_x_ObsMetricValsDict,
        setUpObsCol
    )

    # Flag for using bootstrap sampling
    doBootstrapSampling = False
    numBootstrapSamples = 2

    # Number of metrics to tune on
    numMetricsToTune = 162

    # L1 regularization coefficient, i.e., penalty on param perturbations in objFnc
    # Increase this value to 0.1 or 0.5 or so if you want to reduce 
    # the size of tuned parameter perturbations.
    reglrCoef = 0.0

    # Use these flags to determine whether or not to create specific plots
    createPlotType = {
        'paramsErrorBarsFig': True,
        'biasesOrderedArrowFig': True,
        'threeDotFig': True,
        'metricsBarChart': True,
        'paramsIncrsBarChart': True,
        'paramsAbsIncrsBarChart': True,
        'paramsTotContrbBarChart': True,
        'biasesVsDiagnosticScatterplot': True,
        'dpMin2PtFig': True,
        'dpMinMatrixScatterFig': True,
        'projectionMatrixFigs': True,
        'biasesVsSensMagScatterplot': True,
        'biasesVsSvdScatterplot': True,
        'paramsCorrArrayFig': True,
        'sensMatrixAndBiasVecFig': True,
        'PcaBiplot': True,
        'PcSensMap': True,
        'vhMatrixFig': True,
    }

    # Metrics are observed quantities that we want a tuned simulation to match.
    #    The first column is the metric name.
    #    The order of metricNames determines the order of rows in sensMatrix.
    # The second column is a vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fitting process.
    #   Use a large weight for global (GLB) metrics.
    # The third column is a vector of normalization values for metrics.  
    #   If a value in the 3rd column is set to -999, then the metric is simply normalized by the observed value.
    #   Otherwise, the value in the 3rd column is itself the normalization value for the metric.  
    metricsNamesWeightsAndNormsCustom = \
        [
            #                         ['RESTOM_GLB', 4.0, 10.],
            #                         ['SWCF_GLB', 16.0e-6, -999],
            #                         ['SWCF_DYCOMS', 4.0e-6, -999],
            #                         ['SWCF_HAWAII', 4.00e-6, -999],
            #                         ['SWCF_VOCAL', 4.00e-6, -999],
            #                         ['SWCF_VOCAL_near', 1.00e-6, -999],
            #                         ['SWCF_LBA', 1.00e-6, -999],
            #                         ['SWCF_WP', 1.00e-6, -999],
            #                         ['SWCF_EP', 1.00e-6, -999],
            #                         ['SWCF_NP', 1.00e-6, -999],
            #                         ['SWCF_SP', 1.00e-6, -999],
            #                         ['SWCF_PA', 1.01, -999],
            #                         ['SWCF_CAF', 1.00, -999],
            #                         ['SWCF_Namibia', 4.00e-6, -999],
            #                         ['SWCF_Namibia_near', 1.00e-6, -999],
            #                         ['LWCF_GLB',1.00e-6, -999],
        ]

    # Split up the list above into metric names and the corresponding weights.
    dfMetricsNamesWeightsAndNormsCustom = \
        pd.DataFrame(metricsNamesWeightsAndNormsCustom,
                     columns=['metricsNamesCustom', 'metricsWeightsCustom', 'metricsNormsCustom'])
    metricsNamesCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsNamesCustom']].to_numpy().astype(str)[:, 0]
    metricsWeightsCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsWeightsCustom']].to_numpy().astype(float)
    metricsNormsCustom = dfMetricsNamesWeightsAndNormsCustom[['metricsNormsCustom']].to_numpy().astype(float)

    # These are some metrics that we want to include
    #      in the metrics bar-chart, 3-dot plot, etc.
    # They must be a subset of metricsNames
    highlightedMetricsToPlot = np.array(['SWCF_6_14', 'SWCF_6_18', 'SWCF_8_13',
                                         'SWCF_3_14', 'SWCF_1_14', 'SWCF_3_6', 'SWCF_1_6'])

    # Parameters are tunable model parameters, e.g. clubb_C8.
    # The float listed below after the parameter name is a factor that is used below for scaling plots.
    #   It is not a weight and doesn't affect optimized values; it just makes the plots more readable.
    # Each parameter is associated with two sensitivity simulations; in one, the parameter is perturbed
    #    up and in the other, it is perturbed down.
    #    The output from each sensitivity simulation is expected to be stored in its own netcdf file.
    #    Each netcdf file contains metric values and parameter values for a single simulation.
    folder_name = 'Regional_files/20241022_1yr_20x20regs/20.0sens1022_'

    paramsNamesScalesAndFilenames = \
        [
            ['clubb_c8', 1.0,
             folder_name + '14_Regional.nc',
             folder_name + '15_Regional.nc'],
            ['clubb_c_invrs_tau_n2', 1.0,
             folder_name + '10_Regional.nc',
             folder_name + '11_Regional.nc'],
            ['clubb_c_invrs_tau_sfc', 1.0,
             folder_name + '6_Regional.nc',
             folder_name + '7_Regional.nc'],
            ['clubb_c_invrs_tau_wpxp_n2_thresh', 1.e3,
             folder_name + '8_Regional.nc',
             folder_name + '9_Regional.nc'],
            ['clubb_c_invrs_tau_n2_wp2', 1.0,
             folder_name + '4_Regional.nc',
             folder_name + '5_Regional.nc'],
        ]

    # Split up the above list into parameter names, scales, and filenames.
    dfparamsNamesScalesAndFilenames = \
        pd.DataFrame(paramsNamesScalesAndFilenames,
                     columns=['paramsNames', 'paramsScales',
                              'sensNcFilenames', 'sensNcFilenamesExt'])
    # 'sensNcFilenamesExt', 'sensNcFilenames'] )
    paramsNames = dfparamsNamesScalesAndFilenames[['paramsNames']].to_numpy().astype(str)[:, 0]
    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    paramsScales = dfparamsNamesScalesAndFilenames[['paramsScales']].to_numpy().astype(float)[:, 0]
    sensNcFilenames = dfparamsNamesScalesAndFilenames[['sensNcFilenames']].to_numpy().astype(str)[:, 0]
    sensNcFilenamesExt = dfparamsNamesScalesAndFilenames[['sensNcFilenamesExt']].to_numpy().astype(str)[:, 0]

    # Below we designate the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array([''])

    prescribedParamsNamesScalesAndValues = \
        [
            # ['clubb_c11b', 1.0, 0.5,
            #  folder_name + 'clubb_c11bm_Regional.nc',
            #  folder_name + 'clubb_c11bp_Regional.nc'],
        ]
    # Split up the above list into parameter names, scales, and filenames.
    dfprescribedParamsNamesScalesAndValues = \
        pd.DataFrame(prescribedParamsNamesScalesAndValues,
                     columns=['prescribedParamsNames',
                              'prescribedParamsScales',
                              'prescribedParamVals',
                              'prescribedSensNcFilenames', 'prescribedSensNcFilenamesExt'
                              ]
                     )
    prescribedParamsNames = dfprescribedParamsNamesScalesAndValues[['prescribedParamsNames']].to_numpy().astype(str)[:,
                            0]
    # Extract scaling factors of parameter values from user-defined list paramsNamesScalesAndFilenames.
    # The scaling is not used for any calculations, but it allows us to avoid plotting very large or small values.
    prescribedParamsScales = dfprescribedParamsNamesScalesAndValues[['prescribedParamsScales']].to_numpy().astype(
        float)[:, 0]
    prescribedParamVals = dfprescribedParamsNamesScalesAndValues[['prescribedParamVals']].to_numpy().astype(float)[:, 0]
    prescribedParamValsRow = prescribedParamVals
    prescribedSensNcFilenames = dfprescribedParamsNamesScalesAndValues[['prescribedSensNcFilenames']].to_numpy().astype(
        str)[:, 0]
    prescribedSensNcFilenamesExt = dfprescribedParamsNamesScalesAndValues[
                                       ['prescribedSensNcFilenamesExt']].to_numpy().astype(str)[:, 0]
    prescribedTransformedParamsNames = np.array([''])

    # Netcdf file containing metric and parameter values from the default simulation
    defaultNcFilename = \
        (
                folder_name + '1_Regional.nc'
        )

    # Metrics from the global simulation that uses the tuner-recommended parameter values
    globTunedNcFilename = \
        (
                folder_name + '69_Regional.nc'
            # defaultNcFilename
        )

    # Comment out if not using 20x20reg files
    varPrefixes = ["SWCF"]
    # varPrefixes = ["SWCF", "LWCF", "PRECT"]
    metricsNamesWeightsAndNorms, metricGlobalValsFromFile \
        = setUp_x_MetricsList(varPrefixes, defaultNcFilename)
    # Split up the list above into metric names and the corresponding weights.
    dfMetricsNamesWeightsAndNorms = \
        pd.DataFrame(metricsNamesWeightsAndNorms, columns=['metricsNames', 'metricsWeights', 'metricsNorms'])
    metricsNames = dfMetricsNamesWeightsAndNorms[['metricsNames']].to_numpy().astype(str)[:, 0]
    metricsWeights = dfMetricsNamesWeightsAndNorms[['metricsWeights']].to_numpy().astype(float)
    # metricsNorms = dfMetricsNamesWeightsAndNorms[['metricsNorms']].to_numpy().astype(float)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    metricGlobalAvgs = np.diag(np.dot(metricsWeights.reshape(-1, len(varPrefixes), order='F').T,
                                      defaultMetricValsCol.reshape(-1, len(varPrefixes), order='F')))

    if not np.isclose(metricGlobalValsFromFile, metricGlobalAvgs).all():
        print("Error: metricGlobalAvgs not equal to metricGlobalValsFromFile")
    print("\nThe following two global values should be close to each other:")
    print("metricGlobalAvgs =", metricGlobalAvgs)
    print("metricGlobalValsFromFile =", metricGlobalValsFromFile)
    if beVerbose:
        print("defaultMetricValsCol printed as array = ")
        # Calculate number of regions in the east-west (X) and north-south (Y) directions
        boxSize = 20
        numXBoxes = np.rint(360 / boxSize).astype(int)  # 18
        numYBoxes = np.rint(180 / boxSize).astype(int)  # 9
        defaultMetricValsReshaped = defaultMetricValsCol.reshape((numYBoxes, numXBoxes))
        # defaultMetricValsRolled = np.roll(defaultMetricValsReshaped, -9, axis=1)
        np.set_printoptions(linewidth=200)
        print(np.around(defaultMetricValsReshaped, 2))
        # print(np.around(defaultMetricValsRolled,2))

    # Read observed values of regional metrics on regular tiled grid into a Python dictionary
    (obsMetricValsDict, obsWeightsDict) = \
        (
            setUp_x_ObsMetricValsDict(varPrefixes, suffix="_[0-9]+_", obsPathAndFilename=folder_name + "20241011_20.0_OBS.nc")
        )

    # Set metricsNorms to be a global average
    obsGlobalAvgObsWeights = np.zeros(len(varPrefixes))
    obsGlobalStdObsWeights = np.zeros(len(varPrefixes))
    obsGlobalAvgCol = np.empty(shape=[0, 1])
    obsGlobalStdCol = np.empty(shape=[0, 1])
    for idx, varPrefix in np.ndenumerate(varPrefixes):
        keysVarPrefix = [key for key in obsWeightsDict.keys() if varPrefix in key]
        # obsWeightsNames = np.array(list(obsWeightsDict.keys()), dtype=str)
        obsWeightsNames = np.array(keysVarPrefix, dtype=str)
        obsWeightsUnnormlzd = setUpObsCol(obsWeightsDict, obsWeightsNames)
        obsWeights = obsWeightsUnnormlzd / np.sum(obsWeightsUnnormlzd)
        # metricsWeights = obsWeights
        # obsWeights = np.vstack([obsWeights] * len(varPrefixes))
        metricsNamesVarPrefix = [key for key in obsMetricValsDict.keys() if varPrefix in key]
        obsMetricValsColVarPrefix = setUpObsCol(obsMetricValsDict, metricsNamesVarPrefix)
        obsGlobalStdObsWeights[idx] = np.std(obsMetricValsColVarPrefix)
        obsGlobalAvgObsWeights[idx] = np.dot(obsWeights.T, obsMetricValsColVarPrefix)
        # For sea-level pressure, the global avg is too large to serve as a representative normalization
        if varPrefix == 'PSL':
            obsGlobalAvgObsWeights[idx] = 1e-3 * obsGlobalAvgObsWeights[idx]
        print(f"obsGlobalAvgObsWeights for {varPrefix} =", obsGlobalAvgObsWeights[idx])
        obsGlobalAvgCol = np.vstack((obsGlobalAvgCol,
                                     obsGlobalAvgObsWeights[idx] * np.ones((len(obsWeights), 1))
                                     ))
        obsGlobalStdCol = np.vstack((obsGlobalStdCol,
                                     obsGlobalStdObsWeights[idx] * np.ones((len(obsWeights), 1))
                                     ))
    # Warning: Using a global average as the constant weight produces little normalized
    #     sensitivity for PSL
    metricsNorms = np.copy(obsGlobalAvgCol)
    # metricsNorms = np.copy(obsGlobalStdCol)

    # obsMetricValsReshaped = obsMetricValsCol.reshape((9,18))
    # biasMat = defaultMetricValsReshaped - obsMetricValsReshaped
    # print("biasMat =")
    # print(np.around(biasMat,2))

    # mse = np.sum(metricsWeights*(defaultMetricValsCol - obsMetricValsCol)**2) \
    #   / np.sum(metricsWeights)
    # rmse = np.sqrt(mse)
    # print("rmse between default and obs =", rmse)

    # The special regions are tacked onto the end of
    #     the usual metrics vectors
    numMetricsNoCustom = len(metricsNames)

    metricsNames = np.append(metricsNames, metricsNamesCustom)
    metricsWeights = np.vstack((metricsWeights, metricsWeightsCustom))
    metricsNorms = np.vstack((metricsNorms, metricsNormsCustom))

    metricsNamesNoprefix = np.char.replace(metricsNames, "SWCF_", "")



    # Observed values of our metrics, from, e.g., CERES-EBAF.
    # These observed metrics will be matched as closely as possible by analyzeSensMatrix.
    # NOTE: PRECT is in the unit of m/s
    (obsMetricValsDictCustom, obsWeightsDictCustom) = \
        (
            setUp_x_ObsMetricValsDict(metricsNamesCustom, suffix="", obsPathAndFilename=folder_name + "20241011_20.0_OBS.nc")
        )

    # For custom regions, make simulated values a numpy float,
    #     like the other metrics
    #obsMetricValsDictCustom = {key: np.float32(value) \
    #                           for key, value in obsMetricValsDictCustom.items()}

    # Add obs of custom metrics to obs dictionary
    obsMetricValsDict.update(obsMetricValsDictCustom)

    # Sanity check: is highlightedMetricsToPlot a subset of metricsNames?
    if np.setdiff1d(highlightedMetricsToPlot, metricsNames).size != 0:
        print("One of the metrics names specified in highlightedMetricsToPlot "
              "does not appear in metricsNames:")
        print(np.setdiff1d(highlightedMetricsToPlot, metricsNames))

    return (numMetricsNoCustom,
            metricsNames, metricsNamesNoprefix,
            varPrefixes,
            highlightedMetricsToPlot, createPlotType,
            metricsWeights, metricsNorms,
            obsMetricValsDict,
            paramsNames, paramsScales,
            transformedParamsNames,
            prescribedParamsNames, prescribedParamsScales,
            prescribedTransformedParamsNames,
            prescribedParamValsRow,
            prescribedSensNcFilenames, prescribedSensNcFilenamesExt,
            sensNcFilenames, sensNcFilenamesExt,
            defaultNcFilename, globTunedNcFilename,
            reglrCoef, doBootstrapSampling, numBootstrapSamples, numMetricsToTune)


def abbreviateParamsNames(paramsNames):
    """
    Abbreviate parameter names so that they fit on plots.
    This is handled manually with the lines of code below.
    """

    paramsAbbrv = np.char.replace(paramsNames, 'clubb_', '')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'c_invrs_tau_', '')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'wpxp_n2', 'n2')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'altitude', 'alt')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'threshold', 'thres')
    paramsAbbrv = np.char.replace(paramsAbbrv, 'thresh', 'thres')

    return paramsAbbrv
