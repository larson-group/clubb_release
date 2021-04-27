"""
    Construct and analyze a sensitivity matrix
    that tells us how user-defined metrics change
    with changes in parameter values.
"""

def main():

    import numpy as np
    import pdb

    # Metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsWeights = np.array([[1.], [1.], [1.]])

    # Parameters are tunable model parameters.
    #    The order of paramsNames must match the order of filenames below.
    paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

    # This the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array([''])

    # Netcdf file containing metric and parameter values from the default simulation
    defaultNcFilename = \
        'default.nc' #'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux.ne30_ne30_GLBmean.nc'

    # Metrics from simulation that use the SVD-recommended parameter values
    # Here, we use default simulation just as a placeholder.
    linSolnNcFilename = \
        'default.nc' #'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux.ne30_ne30_GLBmean.nc'

    # This is a list of one netcdf file per each sensitivity simulation.
    # Each file contains metrics and parameter values for a single simulation.
    # There should be one sensitivity simulation per each tunable parameter.
    # These filenames must be listed in the same order as the parameters (paramsNames).
    sensNcFilenames = \
    np.array([ \
        #'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux_c82.ne30_ne30_GLBmean.nc', \
        #'/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux_n21.ne30_ne30_GLBmean.nc' \
        'sens0.nc', \
        'sens1.nc'    
              ])

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    # These observed metrics will be matched as closely as possible by analyzeSensMatrix.
    # Global mean, NOTE PRECT is in the unit of m/s
    obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, \
    defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    defaultParamValsOrigRow, dparamsSoln, paramsSoln, paramsSolnPC = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                      metricsWeights,
                      sensNcFilenames, defaultNcFilename,
                      obsMetricValsDict)

    # See if new global simulation output based on a linear combination
    #    of the SVD-calculated parameter values matches what we expect.
    linSolnBias = calcLinSolnBias(linSolnNcFilename, defaultNcFilename,
                                  metricsNames)

    # Create a heatmap plot that allows us to visualize the normalized sensitivity matrix
    plotNormlzdSensMatrix(normlzdSensMatrix, metricsNames, paramsNames)

    print("\nReached the end of main.")

    return

def analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                      metricsWeights,
                      sensNcFilenames, defaultNcFilename,
                      obsMetricValsDict):
    """
    Input: Information about metrics and parameter values from
            a default simulation, some sensitivity simulations, and observations.
    Output: The inverse of the singular value decomposition.
    """

    import numpy as np
    import sys
    import netCDF4
    import pdb

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    print("\ndefaultNcFilename = ")
    print(defaultNcFilename)

    print("\nsensNcFilenames = ")
    print(sensNcFilenames)

    # Number of tunable parameters
    numParams = len(paramsNames)

    #print("\nnumParams =")
    #print(numParams)

    # Number of metrics
    numMetrics = len(metricsNames)

    #print("\nnumMetrics =")
    #print(numMetrics)

    print("\nmetricsNames = ")
    print(metricsNames)

    print("\nmetricsWeights = ")
    print(metricsWeights)

    print("\nparamsNames = ")
    print(paramsNames)

    print("\ntransformedParamsNames = ")
    print(transformedParamsNames)

    # Set up a column vector of observed metrics
    obsMetricValsCol = \
            setupObsCol(obsMetricValsDict, metricsNames, numMetrics)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow = \
            setupDefaultVectors(metricsNames, paramsNames, transformedParamsNames,
                                numMetrics, numParams,
                                defaultNcFilename)

    # Make sure that each sensitivity simulation changes one and only one parameter.
    # However, here we only check parameters that are listed in paramsNames.
    for sensFileIdx in np.arange(numParams):
        # Read netcdf file of a sensitivity simulation.
        f_sensParams = netCDF4.Dataset(sensNcFilenames[sensFileIdx], 'r')
        # Now loop over parameters and compare values in the sensitivity and default simulations
        for paramIdx in np.arange(numParams):
            paramName = paramsNames[paramIdx]
            # Assume each metric is stored as length-1 array, rather than scalar.
            #   Hence the "[0]" at the end is needed.
            sensParamVal = np.ndarray.item( f_sensParams.variables[paramName][0] )
            # Transform [0,1] variable to extend over range [0,infinity]
            if paramName in transformedParamsNames:
                #sensParamVal = -np.log(1-sensParamVal)
                sensParamVal = np.log(sensParamVal)
            defaultParamVal = defaultParamValsRow[0][paramIdx]
            sensDefaultAreClose = np.isclose(sensParamVal,defaultParamVal)
            if sensFileIdx==paramIdx and sensDefaultAreClose:
                print("\nsensFileIdx =", sensFileIdx)
                print("paramIdx =", paramIdx, ", defaultParamVal =", defaultParamVal, ", sensParamVal =", sensParamVal)
                sys.exit("Error: a sensitivity simulation has left a designated parameter at its default value.")
            if sensFileIdx!=paramIdx and not sensDefaultAreClose:
                print("\nsensFileIdx =", sensFileIdx)
                print("paramIdx =", paramIdx, ", defaultParamVal =", defaultParamVal, ", sensParamVal =", sensParamVal)
                sys.exit("Error: sensitivity simulation has changed the value of an undesignated parameter.")
        f_sensParams.close()

    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow = \
            setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                    numMetrics, numParams,
                    sensNcFilenames)

    # Calculate the magnitude of the maximum value of parameters
    #    from the default and sensitivity runs, for later use
    #    in scaling the normalized sensitivity matrix.
    maxMagParamValsRow = np.maximum( np.abs(defaultParamValsRow), \
                                     np.abs(sensParamValsRow) )
    if np.any( np.isclose(maxMagParamValsRow, np.zeros((1,numParams))) ):
        print("\nmaxMagParamValsRow =")
        print(maxMagParamValsRow)
        sys.exit("Error: A parameter value from both default and sensitivity simulation is zero.")    

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    # Use untransformed (original) parameter values.
    defaultBiasesCol, sensMatrixOrig, normlzdSensMatrixOrig = \
         constructSensMatrix(sensMetricValsMatrix, sensParamValsOrigRow,
                            defaultMetricValsCol, defaultParamValsOrigRow,
                            np.full_like(maxMagParamValsRow,1.0),
                            obsMetricValsCol,
                            numMetrics, numParams,
                            beVerbose=False)

    print("\nsensMatrixOrig =")
    print(sensMatrixOrig)

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    # Use transformed parameter values.
    defaultBiasesCol, sensMatrix, normlzdSensMatrix = \
         constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                            defaultMetricValsCol, defaultParamValsRow,
                            maxMagParamValsRow,
                            obsMetricValsCol,
                            numMetrics, numParams,
                            beVerbose=True)

    # In order to de-weight certain metrics, multiply each row of normlzdSensMatrix 
    # by metricsWeights
    normlzdWeightedSensMatrix = np.diag(np.transpose(metricsWeights)[0]) @ normlzdSensMatrix


    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    svdInvrsNormlzdWeighted, svdInvrsNormlzdWeightedPC, sValsTruncInvNormlzdWeighted, \
    vhNormlzdWeighted = \
    calcSvdInvrs(normlzdWeightedSensMatrix)

    #print("\nNormalized, weighted SVD inverse =")
    #print(svdInvrsNormlzdWeighted)

    paramsSoln, dparamsSoln, defaultBiasesApprox = \
    calcParamsSoln(svdInvrsNormlzdWeighted, metricsWeights, maxMagParamValsRow, \
                   sensMatrix, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeighted, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames )

    paramsSolnPC, dparamsSolnPC, defaultBiasesApproxPC = \
    calcParamsSoln(svdInvrsNormlzdWeightedPC, metricsWeights, maxMagParamValsRow, \
                   sensMatrix, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeighted, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames )

    print("\ndefaultBiasesApprox =")
    print(defaultBiasesApprox)

    print("\ndparamsSoln =")
    print(dparamsSoln)

    print("\nparamsSoln =")
    print(paramsSoln)

    # This check is currently broken if any params span [0,1], e.g., C5
    #if (transformedParamsNames == np.array([''])).all():
    print("\nCheck: Does defaultBiasesOrigApprox approximate defaultBiasesCol above?")
    defaultBiasesOrigApprox = sensMatrixOrig @ dparamsSoln
    print("defaultBiasesOrigApprox =")
    print(defaultBiasesOrigApprox)

    defaultBiasesOrigApproxPC = sensMatrixOrig @ dparamsSolnPC

    return (defaultMetricValsCol, defaultBiasesCol, \
            defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
            sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
            defaultParamValsOrigRow, dparamsSoln, paramsSoln, paramsSolnPC)

def constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                        defaultMetricValsCol, defaultParamValsRow,
                        maxMagParamValsRow,
                        obsMetricValsCol,
                        numMetrics, numParams,
                        beVerbose):
    """
    Inputs: Metric and parameter values from default and sensitivity simulations,
        and from observations.
    Output: A sensitivity matrix listing the partial derivatives dmetric/dparam.
    """

    import numpy as np
    import sys
    import pdb

    # Matrix of metric values from default simulation
    # Each column in the matrix is repeated numParams times, for later multiplication
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    #print("\ndefaultMetricValsMatrix =")
    #print(defaultMetricValsMatrix)

    # Sensitivity simulation metrics minus default simulation metrics
    dmetricsMatrix = np.subtract(sensMetricValsMatrix, defaultMetricValsMatrix)

    if beVerbose:
        print("\ndmetricsMatrix =")
        print(dmetricsMatrix)

    # Perturbation parameter values
    dparamsRow = np.subtract(sensParamValsRow, defaultParamValsRow)

    if beVerbose:
        print("\ndparamsRow =")
        print(dparamsRow)

    # Make sure that the parameter values from sensitivity simulations
    #    are actually different than the ones from the default simulation:
    if np.any( np.isclose(dparamsRow, np.zeros((1,numParams))) ):
        print("\ndparamsRow =")
        print(dparamsRow)
        sys.exit("Error: A sensitivity simulation has left its 'changed' parameter at the default value.")

    #print("\nreciprocal of dparamsRow")
    #print(np.reciprocal(dparamsRow))

    # Matrix of inverse perturbation parameter values.
    # Used for forming sensitivity derivatives.
    invrsDparamsMatrix = np.ones((numMetrics,1)) @ np.reciprocal(dparamsRow)

    #print("\ninvrsDparamsMatrix =")
    #print(invrsDparamsMatrix)

    # Sensitivity matrix of derivatives, dmetrics/dparams
    sensMatrix = dmetricsMatrix * invrsDparamsMatrix

    if beVerbose:
        print("\nsensMatrix =")
        print(sensMatrix)

    # Store biases in default simulation
    defaultBiasesCol = np.subtract(obsMetricValsCol, defaultMetricValsCol)

    if beVerbose:
        print("\ndefaultBiasesCol =")
        print(defaultBiasesCol)

    # Matrix of inverse biases.
    # Used for forming normalized sensitivity derivatives.
    invrsBiasesMatrix = np.reciprocal(defaultBiasesCol) @ np.ones((1,numParams))

    #print("\ninvrsBiasesMatrix =")
    #print(invrsBiasesMatrix)

    # Form matrix of default parameter values, for later normalization of the sensitivity matrix
    maxMagParamValsMatrix = np.ones((numMetrics,1)) @ maxMagParamValsRow

    # Sensitivity matrix, normalized by biases and parameter values
    normlzdSensMatrix = sensMatrix * invrsBiasesMatrix * maxMagParamValsMatrix

    if beVerbose:
        print("\nnormlzdSensMatrix =")
        print(normlzdSensMatrix)

    return  (defaultBiasesCol, sensMatrix, normlzdSensMatrix)


def calcSvdInvrs(normlzdWeightedSensMatrix):
    """
    Input: sensitivity matrix
    Output: singular value decomposition of sensitivity matrix
        and related quantities, like the inverse of the SVD.
    """


    import numpy as np
    import sys
    import pdb

    u, s, vh = np.linalg.svd(normlzdWeightedSensMatrix, full_matrices=False)

    print("\nSingular values =")
    print(s)

    #print("\nvh =")
    #print(vh)

    #print("\nu =")
    #print(u)

    sValsTrunc = s

    sValsTrunc[sValsTrunc < 1e-10] = -1

    print("\nTruncated Singular values =")
    print(sValsTrunc)

    sValsTruncInv = np.reciprocal(sValsTrunc)

    sValsTruncInv[sValsTruncInv < 0] = 0

    #print("\nInverse truncated singular values =")
    #print(sValsTruncInv)

    svdInvrs = np.transpose(vh) @ np.diag(sValsTruncInv) @ np.transpose(u)

    # Assertion check: is svdInvrs truly the inverse of normlzdSensMatrix?
    numParams = normlzdWeightedSensMatrix.shape[1] # = number of columns
    if not np.all( np.isclose(np.identity(numParams), svdInvrs @ normlzdWeightedSensMatrix, \
                      rtol=1e-6 , atol=1e-6 ) ):
        sys.exit("Error: svdInvrs is not the inverse of normlzdWeightedSensMatrix")

    #print("\nSVD inverse =")
    #print(svdInvrs)


    # Delete half the singular values in order to show just the most important patterns.
    #pdb.set_trace()
    halfNumSVals = np.floor_divide( sValsTrunc.size, 2 )
    sValsInvPC = np.zeros_like(sValsTruncInv)
    sValsInvPC[0:halfNumSVals] = sValsTruncInv[0:halfNumSVals]
    #sValsInvPC[0:sValsTruncInv.size-3] = sValsTruncInv[0:sValsTruncInv.size-3]

    #pdb.set_trace()

    svdInvrsPC = np.transpose(vh) @ np.diag(sValsInvPC) @ np.transpose(u)

    eigVals, eigVecs = np.linalg.eig(np.transpose(normlzdWeightedSensMatrix) @ normlzdWeightedSensMatrix)

    #print("\neigVals =")
    #print(eigVals)

    #print("\neigVecs = ")
    #print(eigVecs)

    return ( svdInvrs, svdInvrsPC, sValsTruncInv, vh )

def calcParamsSoln(svdInvrsNormlzdWeighted, metricsWeights, maxMagParamValsRow, \
                   sensMatrix, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeighted, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames ):

    import numpy as np
    import pdb

    # Calculate solution in transformed space
    dparamsSoln = svdInvrsNormlzdWeighted @ metricsWeights * np.transpose(maxMagParamValsRow)
    defaultBiasesApprox = sensMatrix @ dparamsSoln
    paramsSoln = np.transpose(defaultParamValsOrigRow) + dparamsSoln
    paramsErrMatrix = np.transpose(vhNormlzdWeighted) * sValsTruncInvNormlzdWeighted
    paramsErrStd = np.sqrt( np.reshape( np.sum( paramsErrMatrix * paramsErrMatrix, axis=1 ), (-1, 1)) )
    #pdb.set_trace()
    # Transform some variables from [-inf,inf] back to [0,inf] range
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        if paramName in transformedParamsNames:
            #paramsSoln[idx,0] = 1.0-np.exp(-paramsSoln[idx,0])
            paramsSoln[idx,0] = np.exp(dparamsSoln[idx,0]) * defaultParamValsOrigRow[0,idx]
            dparamsSoln[idx,0] = paramsSoln[idx,0] - defaultParamValsOrigRow[0,idx]

    return ( paramsSoln, dparamsSoln, defaultBiasesApprox )

def setupObsCol(obsMetricValsDict, metricsNames, numMetrics):
    """
    Input: A python dictionary of observed metrics.
    Output: A column vector of observed metrics
    """

    import numpy as np
    import pdb

    # Set up column vector of numMetrics elements containing
    # "true" metric values from observations
    obsMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        obsMetricValsCol[idx] = obsMetricValsDict[metricName]

    print("\nobsMetricValsCol =")
    print(obsMetricValsCol)

    return obsMetricValsCol

def setupDefaultVectors(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
                        defaultNcFilename):
    """
    Input: Filename containing default-simulation metrics and parameters.
    Output: Column vector of default-sim metrics,
            and row vector of default-sim parameter values.
    """

    import numpy as np
    import netCDF4

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

    print("\ndefaultMetricValsCol =")
    print(defaultMetricValsCol)

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

    print("\ndefaultParamValsOrigRow =")
    print(defaultParamValsOrigRow)

    print("\ndefaultParamValsRow =")
    print(defaultParamValsRow)

    f_defaultMetricsParams.close()

    return (defaultMetricValsCol, defaultParamValsRow, defaultParamValsOrigRow)

def setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                    numMetrics, numParams,
                    sensNcFilenames):
    """
    Input: List of filenames, one per each sensitivity simulation.
    Output: Row vector of modified parameter values from sensitivity simulations.
            Matrix of metrics, where each column corresponds to
                a single sensitivity simulation.
    """


    import numpy as np
    import netCDF4

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

    #sensParamValsRow = np.array([[2., 4.]])

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
            sensMetricValsMatrix[row,col] = f_sens.variables[metricName][0]
        f_sens.close()

    #sensMetricValsMatrix = np.array([[1., 2.], [3., 4.]])

    print("\nsensMetricValsMatrix =")
    print(sensMetricValsMatrix)

    return(sensMetricValsMatrix, sensParamValsRow, sensParamValsOrigRow)

def plotNormlzdSensMatrix(normlzdSensMatrix, metricsNames, paramsNames):

    import numpy as np
    import matplotlib.pyplot as plt    

    fig, ax = plt.subplots(1,1)
    coloredMatrix = ax.imshow(normlzdSensMatrix, aspect='auto', cmap=plt.cm.RdBu)
    # Print the matrix element values as text
    for (i, j), value in np.ndenumerate(normlzdSensMatrix):
        plt.text(j, i, "%.3f"%value, va='center', ha='center')
    ax.set_yticks(np.arange(metricsNames.size))
    ax.set_yticklabels(metricsNames)
    ax.set_xticks(np.arange(paramsNames.size))
    ax.set_xticklabels(paramsNames)
    plt.ion()  # Prevent plotting from blocking rest of execution
    plt.show()
    plt.pause(1)

def calcLinSolnBias(linSolnNcFilename, defaultNcFilename,
                    metricsNames):

    import numpy as np
    import netCDF4

    # Number of metrics
    numMetrics = len(metricsNames)

    # Read netcdf file with metrics and parameters from default simulation
    f_default = netCDF4.Dataset(defaultNcFilename, 'r')
    # Set up column vector of numMetrics elements containing
    # metric values from default simulation
    defaultMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        defaultMetricValsCol[idx] = f_default.variables[metricName][0]
    f_default.close()

    # Read netcdf file with metrics and parameters from default simulation
    f_linSoln = netCDF4.Dataset(linSolnNcFilename, 'r')
    # Set up column vector of numMetrics elements containing
    # metric values from simulation that uses SVD-recommended parameter values
    linSolnMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        linSolnMetricValsCol[idx] = f_linSoln.variables[metricName][0]
    f_linSoln.close()

    linSolnBias = linSolnMetricValsCol - defaultMetricValsCol

    print("\nlinSolnBias")
    print(linSolnBias)

    return linSolnBias

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
