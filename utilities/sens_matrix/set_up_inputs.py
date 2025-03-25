# -*- coding: utf-8 -*-

# Run this app with `python3 quadtune_driver.py` and
# view the plots at http://127.0.0.1:8050/ in your web browser.
# (To open a web browser on a larson-group computer,
# login to malan with `ssh -X` and then type `firefox &`.)

"""
This file contains a set of utility functions that read
the netcdf files containing regional metric values, observations,
and parameter values, and set up various python arrays,
including the default bias column vector and the parameter
value row vector.
"""

import numpy as np
from re import search, match
import netCDF4
import sys



def setUpColAndRowVectors(metricsNames, metricsNorms,
                          obsMetricValsDict,
                          paramsNames, transformedParamsNames,
                          prescribedParamsNames, prescribedParamValsRow,
                          prescribedTransformedParamsNames,
                          sensNcFilenames,
                          defaultNcFilename
                          ):
    """
    Given netcdf files that contain output from the default and sensitivity global simulations,
    set up column vectors of observations and default-simulation biases.  Also set up
    a row vector of default-simulation parameter values.
    """

    # Set up a column vector of observed metrics
    obsMetricValsCol = setUpObsCol(obsMetricValsDict, metricsNames)

    # Set up a normalization vector for metrics, normMetricValsCol.
    # It equals the observed value when metricsNorms has the special value of -999, 
    #     but otherwise it is set manually in metricsNorms itself.
    normMetricValsCol = np.copy(metricsNorms)
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
    #                  =  default simulation - observations
    defaultBiasesCol = np.subtract(defaultMetricValsCol, obsMetricValsCol)

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

    dnormlzdPrescribedParams = dnormlzdPrescribedParams.T


    return ( obsMetricValsCol, normMetricValsCol,
             defaultBiasesCol,
             defaultParamValsOrigRow,
             magParamValsRow,
             dnormlzdPrescribedParams,
             magPrescribedParamValsRow
           )


def setUp_x_ObsMetricValsDict(varPrefixes, suffix="", obsPathAndFilename=""):
    """
    This is intended for the case in which the metrics are tiles that cover the globe, not custom regions.
    Input: Filename containing observed values of metrics.
    Output: Dictionary of observations.
    """

    # Read netcdf file with metrics and parameters from default simulation
    f_obs = netCDF4.Dataset(obsPathAndFilename, 'r')

    obsMetricValsDict = {}
    obsWeightsDict = {}

    #varPrefixes = ["SWCF"]
    for varName in f_obs.variables:
        #print(varName)
        #         or re.search("^LWCF_[0-9]+_",varName):
        for varPrefix in varPrefixes:
            #if search(f"^{varPrefix}_[0-9]+_", varName):
            if search(f"^{varPrefix}{suffix}",varName):
                #and not "MSWCF" in varName
                varEntry = f_obs[varName]
                varVal = varEntry[:].data[:][0]
                obsMetricValsDict[varName] = varVal
                #print((varName, varVal))
            # Extract observational weights,
            #     which are effectively numpy scalars (0d arrays)
            if search(f"^weights_[0-9]+_[0-9]+_{varPrefix}",varName):
                weightsEntry = f_obs[varName]
                weightsVal = weightsEntry[:].data
                obsWeightsDict[varName] = weightsVal

    f_obs.close()

    #print(obsMetricValsDict)
    #print(len(obsMetricValsDict))

    return (obsMetricValsDict, obsWeightsDict)


def setUp_x_MetricsList(varPrefixes, defPathAndFilename):
    """
    This is intended for the case in which the metrics are tiles that cover the globe, not custom regions.
    Input: Filename of default simulation.
    Output: List of 20x20reg metric values.
    """

    # Read netcdf file with metrics and parameters from default simulation
    f_def = netCDF4.Dataset(defPathAndFilename, 'r')

    metricsNamesWeightsAndNorms = []
    for varPrefix in varPrefixes:
        for varName in f_def.variables:
            #print(varName)
            if match("^numb_[0-9]+_[0-9]+",varName):
                areaWeightEntry = f_def[varName]
                areaWeightVal = areaWeightEntry[:].data[:][0]
                varFullString = varName.replace("numb", varPrefix)
                metricsNamesWeightsAndNorms.append([varFullString,  areaWeightVal, -999])
                #print((SWCF_string, areaWeightVal))

    metricGlobalValsFromFile = np.zeros(len(varPrefixes))
    for index, varPrefix in np.ndenumerate(varPrefixes):
        metricGlobalName = varPrefix + "_GLB"
        metricGlobalValsFromFile[index] = f_def.variables[metricGlobalName][0]

    f_def.close()

    #print(obsMetricValsDict)
    #print(metricsNamesWeightsAndNorms)

    return (metricsNamesWeightsAndNorms, metricGlobalValsFromFile)


def setUpObsCol(obsMetricValsDict, metricsNames):
    """ 
    Input: A python dictionary of observed metrics.
    Output: A column vector of observed metrics
    """

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
    Output: Row vector of default-simulation parameter values.
    """

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


def setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                    numMetrics, numParams,
                    sensNcFilenames,
                    beVerbose):
    """
    Input: List of filenames, one per each sensitivity simulation.
    Output: Row vector of modified parameter values from sensitivity simulations.
            Sensitivity matrix of regional metrics to parameter values, where each column corresponds to
                a single sensitivity simulation.
    """

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
        sensParamValsOrigRow[0, idx] = f_sensParams.variables[paramName][0]
        # Transform [0,1] variable to extend over range [0,infinity]
        if paramName in transformedParamsNames:
            # sensParamValsRow[0,idx] = -np.log(1-sensParamValsRow[0,idx])
            sensParamValsRow[0, idx] = np.log(sensParamValsOrigRow[0, idx])
        else:
            sensParamValsRow[0, idx] = sensParamValsOrigRow[0, idx]
        f_sensParams.close()

    # sensParamValsRow = np.array([[2., 4.]])

    if beVerbose:
        print("\nsensParamValsOrigRow =")
        print(sensParamValsOrigRow)
        print("\nsensParamValsRow =")
        print(sensParamValsRow)

    # numMetrics x numParams matrix of metric values
    # from sensitivity simulations
    sensMetricValsMatrix = np.zeros((numMetrics, numParams))
    for col in np.arange(numParams):
        f_sens = netCDF4.Dataset(sensNcFilenames[col], 'r')
        for row in np.arange(numMetrics):
            metricName = metricsNames[row]
            sensMetricValsMatrix[row, col] = f_sens.variables[metricName][0]
            # if (metricName[0:3] == 'PSL'):  # subtract 9e4 from sea-level pressure for better scaling
            #    sensMetricValsMatrix[row,col] = f_sens.variables[metricName][0] - 9.e4
            #    print("metricName[0:3]=", metricName[0:3])
            #    print("sensMetricValsMatrix=", sensMetricValsMatrix[row][col])
            # else:
            #    sensMetricValsMatrix[row,col] = f_sens.variables[metricName][0]
        f_sens.close()

    # sensMetricValsMatrix = np.array([[1., 2.], [3., 4.]])

    if beVerbose:
        print("\nsensMetricValsMatrix =")
        print(sensMetricValsMatrix)

    return (sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow)


def setupDefaultMetricValsCol(metricsNames, defaultNcFilename):
    """
    Input: Filename containing default-simulation metrics.
    Output: Column vector of default-sim metrics.
    """

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

