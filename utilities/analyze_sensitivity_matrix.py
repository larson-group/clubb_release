"""
    Construct and analyze a sensitivity matrix
    that tells us how user-defined metrics change
    with changes in parameter values.
"""

def main():

    import numpy as np
    import pdb

    # Metrics are observed quantities that we want to match.
    metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

    # Parameters are tunable model parameters.
    #paramsNames = np.array(['C5','C8'])
    paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

    transformedParams = np.array([''])

    # This is a list of one netcdf file per each sensitivity simulation.
    # Each file contains metrics and parameter values for a single simulation.
    # There should be one sensitivity simulation per each tunable parameter.
    # These filenames must be listed in the same order as the parameters (paramsNames).
    #sensNcFilenames = np.array(['/home/vlarson/canopy/scripts/sens_C5_MetricsParams.nc', \
    #                            '/home/vlarson/canopy/scripts/sens_C8_MetricsParams.nc'])
    sensNcFilenames = \
    np.array([ \
        '/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux_c82.ne30_ne30_GLBmean.nc', \
        '/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux_n21.ne30_ne30_GLBmean.nc' \
             ])

    # Netcdf file containing metric and parameter values from the default simulation
    #defaultNcFilename = '/home/vlarson/canopy/scripts/defaultMetricsParams.nc'
    defaultNcFilename = \
        '/home/vlarson/canopy/scripts/anvil.c689c7e.repeatbmg_flux.ne30_ne30_GLBmean.nc'

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}    # Global mean, NOTE PRECT is in the unit of m/s
    #obsMetricValsDict = {'LWCF': 17.39229, 'PRECT': 5.540433451401774e-09, 'SWCF': -63.05767 }    # SC(GPCI) mean, NOTE PRECT is in the unit of m/s

    # Calculate changes in parameter values needed to match metrics.
    sensMatrix, normlzdSensMatrix, svdInvrsNormlzd, dparamsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParams,
                      sensNcFilenames, defaultNcFilename,
                      obsMetricValsDict)

    print("\nReached the end of main.")

    return

def analyzeSensMatrix(metricsNames, paramsNames, transformedParams,
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
    
    print("\nparamsNames = ")
    print(paramsNames)    
    
    # Set up a column vector of observed metrics
    obsMetricValsCol = \
            setupObsCol(obsMetricValsDict, metricsNames, numMetrics)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultMetricValsCol, defaultParamValsRow = \
            setupDefaultVectors(metricsNames, paramsNames, transformedParams,
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
            if paramName in transformedParams:
                sensParamVal = -np.log(1-sensParamVal)
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
    sensParamValsRow, sensMetricValsMatrix = \
            setupSensArrays(metricsNames, paramsNames, transformedParams,
                    numMetrics, numParams,
                    sensNcFilenames)

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    sensMatrix, normlzdSensMatrix = \
         constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                            defaultMetricValsCol, defaultParamValsRow,
                            obsMetricValsCol,
                            numMetrics, numParams)

    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    svdInvrsNormlzd = calcSvdInvrs(normlzdSensMatrix)

    print("\nNormalized SVD inverse =")
    print(svdInvrsNormlzd)

    dparamsSoln = svdInvrsNormlzd @ np.ones((numMetrics,1)) * np.transpose(defaultParamValsRow)
    paramsSoln = np.transpose(defaultParamValsRow) + dparamsSoln
    # Transform some variables from [0,infinity] back to [0,1] range
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        if paramName in transformedParams:
            paramsSoln[idx,0] = 1.0-np.exp(-paramsSoln[idx,0])
            dparamsSoln[idx,0] = 1.0-np.exp(-dparamsSoln[idx,0])

    print("\ndparamsSoln =")
    print(dparamsSoln)

    print("\nparamsSoln =")
    print(paramsSoln)

    return (sensMatrix, normlzdSensMatrix, svdInvrsNormlzd, dparamsSoln)

def constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                        defaultMetricValsCol, defaultParamValsRow,
                        obsMetricValsCol,
                        numMetrics, numParams):
    """
    Inputs: Metric and parameter values from default and sensitivity simulations,
        and from observations.
    Output: A sensitivity matrix listing the partial derivatives dmetric/dparam.
    """

    import numpy as np
    import sys

    # Matrix of metric values from default simulation
    # Each column in the matrix is repeated numParams times, for later multiplication
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    print("\ndefaultMetricValsMatrix =")
    print(defaultMetricValsMatrix)

    # Sensitivity simulation metrics minus default simulation metrics
    dmetricsMatrix = np.subtract(sensMetricValsMatrix, defaultMetricValsMatrix)

    print("\ndmetricsMatrix =")
    print(dmetricsMatrix)

    # Perturbation parameter values
    dparamsRow = np.subtract(sensParamValsRow, defaultParamValsRow)

    print("\ndparamsRow =")
    print(dparamsRow)

    # Make sure that the parameter values from sensitivity simulations
    #    are actually different than the ones from the default simulation:
    if np.any( np.isclose(dparamsRow, np.zeros((1,numParams))) ):
        print("\ndparamsRow =")
        print(dparamsRow)
        sys.exit("Error: A sensitivity simulation has left its 'changed' parameter at the default value.")

    print("\nreciprocal of dparamsRow")
    print(np.reciprocal(dparamsRow))

    # Matrix of inverse perturbation parameter values.
    # Used for forming sensitivity derivatives.
    invrsDparamsMatrix = np.ones((numMetrics,1)) @ np.reciprocal(dparamsRow)

    print("\ninvrsDparamsMatrix =")
    print(invrsDparamsMatrix)


    # Sensitivity matrix of derivatives, dmetrics/dparams
    sensMatrix = dmetricsMatrix * invrsDparamsMatrix

    print("\nsensMatrix =")
    print(sensMatrix)


    # Store biases in default simulation
    defaultBiasesCol = np.subtract(obsMetricValsCol, defaultMetricValsCol)

    print("\ndefaultBiasesCol =")
    print(defaultBiasesCol)

    # Matrix of inverse biases.
    # Used for forming normalized sensitivity derivatives.
    invrsBiasesMatrix = np.reciprocal(defaultBiasesCol) @ np.ones((1,numParams))

    print("\ninvrsBiasesMatrix =")
    print(invrsBiasesMatrix)

    # Form matrix of default parameter values, for later normalization of the sensitivity matrix
    defaultParamValsMatrix = np.ones((numMetrics,1)) @ defaultParamValsRow

    # Sensitivity matrix, normalized by biases and parameter values
    normlzdSensMatrix = sensMatrix * invrsBiasesMatrix * defaultParamValsMatrix

    print("\nnormlzdSensMatrix =")
    print(normlzdSensMatrix)

    return  (sensMatrix, normlzdSensMatrix)


def calcSvdInvrs(normlzdSensMatrix):
    """
    Input: sensitivity matrix
    Output: singular value decomposition of sensitivity matrix
        and related quantities, like the inverse of the SVD.
    """


    import  numpy as np

    u, s, vh = np.linalg.svd(normlzdSensMatrix, full_matrices=False)

    print("\nSingular values =")
    print(s)

    print("\nvh =")
    print(vh)

    print("\nu =")
    print(u)

    sValsTrunc = s

    sValsTrunc[sValsTrunc < 1e-10] = -1

    print("\nTruncated Singular values =")
    print(sValsTrunc)

    sValsTruncInv = np.reciprocal(sValsTrunc)

    sValsTruncInv[sValsTruncInv < 0] = 0

    print("\nInverse truncated singular values =")
    print(sValsTruncInv)

    svdInvrs = np.transpose(vh) @ np.diag(sValsTruncInv) @ np.transpose(u)

    #print("\nSVD inverse =")
    #print(svdInvrs)

    eigVals, eigVecs = np.linalg.eig(np.transpose(normlzdSensMatrix) @ normlzdSensMatrix)

    #print("\neigVals =")
    #print(eigVals)

    #print("\neigVecs = ")
    #print(eigVecs)

    return svdInvrs

def setupObsCol(obsMetricValsDict, metricsNames, numMetrics):
    """
    Input: A python dictionary of observed metrics.
    Output: A column vector of observed metrics
    """

    import numpy as np

    # Set up column vector of numMetrics elements containing
    # "true" metric values from observations
    obsMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        obsMetricValsCol[idx] = obsMetricValsDict[metricName]

    print("\nobsMetricValsCol =")
    print(obsMetricValsCol)

    return obsMetricValsCol

def setupDefaultVectors(metricsNames, paramsNames, transformedParams,
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
    defaultParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        defaultParamValsRow[0,idx] = f_defaultMetricsParams.variables[paramName][0]
        # Transform [0,1] variable to extend over range [0,infinity]
        if paramName in transformedParams:
            defaultParamValsRow[0,idx] = -np.log(1-defaultParamValsRow[0,idx])

    print("\ndefaultParamValsRow =")
    print(defaultParamValsRow)

    f_defaultMetricsParams.close()

    return (defaultMetricValsCol, defaultParamValsRow)

def setupSensArrays(metricsNames, paramsNames, transformedParams,
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
    sensParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        # Read netcdf file with changed parameter values from all sensitivity simulations.
        f_sensParams = netCDF4.Dataset(sensNcFilenames[idx], 'r')
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        sensParamValsRow[0,idx] = f_sensParams.variables[paramName][0]
        # Transform [0,1] variable to extend over range [0,infinity]
        if paramName in transformedParams:
            sensParamValsRow[0,idx] = -np.log(1-sensParamValsRow[0,idx])
        f_sensParams.close()

    #sensParamValsRow = np.array([[2., 4.]])

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

    return(sensParamValsRow, sensMetricValsMatrix)

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
    plt.show()
    
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
