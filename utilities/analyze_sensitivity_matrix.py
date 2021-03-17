"""
    Construct and analyze a sensitivity matrix
    that tells us how user-defined metrics change
    with changes in parameter values.
"""

def main():

    import numpy as np
    import pdb

    # Metrics are observed quantities that we want to match.
    metricsNames = np.array(['SWCF', 'LWCF'])

    # Parameters are tunable model parameters.
    paramsNames = np.array(['C5','C8'])
    #paramsNames = np.array(['clubb_c_invrs_tau_n2','clubb_c_invrs_tau_n2_wp2'])

    # This is a list of one netcdf file per each sensitivity simulation.
    # Each file contains metrics and parameter values for a single simulation.
    # There should be one sensitivity simulation per each tunable parameter.
    # These filenames must be listed in the same order as the parameters (paramsNames).
    sensNcFilenames = np.array(['/home/vlarson/canopy/scripts/sens_C5_MetricsParams.nc', \
                                '/home/vlarson/canopy/scripts/sens_C8_MetricsParams.nc'])
    #sensNcFilenames = \
    #np.array(['/home/vlarson/canopy/scripts/anvil-centos7.manu.repeat.ne30_ne30_GLBmean.nc', \
    #                            '/home/vlarson/canopy/scripts/anvil-centos7.manu.repeat.ne30_ne30_GLBmean.nc'])

    # Netcdf file containing metric and parameter values from the default simulation
    defaultNcFilename = '/home/vlarson/canopy/scripts/defaultMetricsParams.nc'
    #defaultNcFilename = '/home/vlarson/canopy/scripts/anvil-centos7.manu.repeat.ne30_ne30_GLBmean.nc'

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    obsMetricValsDict = {'LWCF': 4., 'PRECT': -999., 'SWCF': 6.}

    # Calculate changes in parameter values needed to match metrics.
    svdInvrs = \
        analyzeSensMatrix(metricsNames, paramsNames,
                      sensNcFilenames, defaultNcFilename,
                      obsMetricValsDict)

    print("\nReached the end.")

    return

def analyzeSensMatrix(metricsNames, paramsNames,
                      sensNcFilenames, defaultNcFilename,
                      obsMetricValsDict):
    """
    Input: Information about metrics and parameter values from
            a default simulation, some sensitivity simulations, and observations.
    Output: The inverse of the singular value decomposition.
    """

    import numpy as np

    if ( len(paramsNames) != len(sensNcFilenames)   ):
        print("Number of parameters must equal number of netcdf files.")
        quit()

    # Values of metrics from the default simulation.  Just for testing.
    defaultMetricValsDict = {'SWCF': 5000., 'LWCF': 6000., 'PRECT': -999000.}

    # Parameter values from default simulation.  Just for testing.
    defaultParamValsDict = {'C15': 42000., 'C5': 1000., 'C8': 2000.}

    # Altered parameter values from sensitivity simulations.  Just for testing.
    sensParamValsDict = {'C5': 2000., 'C15': 42000., 'C8': 4000.}

    # Matrix of metrics from sensitivity simulations.
    # Each column corresponds to one simulation.  Just for testing.
    sensMetricValsRefMatrix = np.array([[1000., 2000.], [3000., 4000.]])

    # Number of tunable parameters
    numParams = len(paramsNames)

    print("\nnumParams =")
    print(numParams)

    # Number of metrics
    numMetrics = len(metricsNames)

    print("\nnumMetrics =")
    print(numMetrics)

    # Set up a column vector of observed metrics
    obsMetricValsCol = \
            setupObsCol(obsMetricValsDict, metricsNames, numMetrics)

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultMetricValsCol, defaultParamValsRow = \
            setupDefaultVectors(defaultMetricValsDict, metricsNames,
                                defaultParamValsDict, paramsNames,
                                numParams, numMetrics,
                                defaultNcFilename)
    # Based on the numParams sensitivity simulations,
    #    set up a row vector of modified parameter values.
    # Also set up numMetrics x numParams matrix,
    #    each column of which lists the metrics
    #    from one of the sensitivity simulations
    sensParamValsRow, sensMetricValsMatrix = \
            setupSensArrays(sensMetricValsRefMatrix, metricsNames,
                    sensParamValsDict, paramsNames,
                    numParams, numMetrics,
                     sensNcFilenames)

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    sensMatrix, normlzdSensMatrix = \
         constructSensMatrix(sensParamValsRow, sensMetricValsMatrix,
                            defaultParamValsRow, defaultMetricValsCol,
                            obsMetricValsCol,
                            numParams, numMetrics)

    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    svdInvrs = calcSvd(normlzdSensMatrix)

    return svdInvrs

def constructSensMatrix(sensParamValsRow, sensMetricValsMatrix,
                        defaultParamValsRow, defaultMetricValsCol,
                        obsMetricValsCol,
                        numParams, numMetrics):
    """
    Inputs: Metric and parameter values from default and sensitivity simulations,
        and from observations.
    Output: A sensitivity matrix listing the partial derivatives dmetric/dparam.
    """

    import numpy as np

    # Matrix of metric values from default simulation
    # Each column in the matrix is repeated numParams times, for later multiplication
    defaultMetricValsMatrix = np.matmul(defaultMetricValsCol, np.ones((1,numParams)))

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

    print("\nreciprocal of dparamsRow")
    print(np.reciprocal(dparamsRow))

    # Matrix of inverse perturbation parameter values.
    # Used for forming sensitivity derivatives.
    invrsDparamsMatrix = np.matmul(np.ones((numMetrics,1)),np.reciprocal(dparamsRow))

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
    invrsBiasesMatrix = np.matmul(np.reciprocal(defaultBiasesCol),np.ones((1,numMetrics)))

    print("\ninvrsBiasesMatrix =")
    print(invrsBiasesMatrix)

    # Sensitivity matrix, normalized by the biases
    normlzdSensMatrix = sensMatrix * invrsBiasesMatrix

    print("\nnormlzdSensMatrix =")
    print(normlzdSensMatrix)

    return  (sensMatrix, normlzdSensMatrix)

def calcSvd(sensMatrix):
    """
    Input: sensitivity matrix
    Output: singular value decomposition of sensitivity matrix
        and related quantities, like the inverse of the SVD.
    """


    import  numpy as np

    u, s, vh = np.linalg.svd(sensMatrix)

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

    svdInvrs = np.transpose(vh)*np.diag(sValsTruncInv) *np.transpose(u)

    print("\nSVD inverse =")
    print(svdInvrs)

    eigVals, eigVecs = np.linalg.eig(np.transpose(sensMatrix)*sensMatrix)

    print("\neigVals =")
    print(eigVals)

    print("\neigVecs = ")
    print(eigVecs)

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

def setupDefaultVectors(defaultMetricValsDict, metricsNames,
                        defaultParamValsDict, paramsNames,
                        numParams, numMetrics,
                        defaultNcFilename):
    """
    Input: Filename containing default-simulation metrics and parameters.
    Output: Column vector of default-sim metrics,
            and row vector of default-sim parameter values.
    """

    import numpy as np
    from scipy.io import netcdf

    # Read netcdf file with metrics and parameters from default simulation
    f_defaultMetricsParams = netcdf.netcdf_file(defaultNcFilename, 'r')

    # Set up column vector of numMetrics elements containing
    # metric values from default simulation
    defaultMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
        metricName = metricsNames[idx]
        #defaultMetricValsCol[idx] = defaultMetricValsDict[metricName]
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
        #defaultParamValsRow[0,idx] = defaultParamValsDict[paramName]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        defaultParamValsRow[0,idx] = f_defaultMetricsParams.variables[paramName][0]

    #defaultParamValsRow = np.array([[1., 2.]])

    print("\ndefaultParamValsRow =")
    print(defaultParamValsRow)

    f_defaultMetricsParams.close()

    return (defaultMetricValsCol, defaultParamValsRow)

def setupSensArrays(sensMetricValsRefMatrix, metricsNames,
                    sensParamValsDict, paramsNames,
                    numParams, numMetrics,
                    sensNcFilenames):
    """
    Input: List of filenames, one per each sensitivity simulation.
    Output: Row vector of modified parameter values from sensitivity simulations.
            Matrix of metrics, where each column corresponds to
                a single sensitivity simulation.
    """


    import numpy as np
    from scipy.io import netcdf

    # Create row vector size numParams containing
    # parameter values from sensitivity simulations
    sensParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        # Read netcdf file with changed parameter values from all sensitivity simulations.
        f_sensParams = netcdf.netcdf_file(sensNcFilenames[idx], 'r')
        #sensParamValsRow[0,idx] = sensParamValsDict[paramName]
        # Assume each metric is stored as length-1 array, rather than scalar.
        #   Hence the "[0]" at the end is needed.
        sensParamValsRow[0,idx] = f_sensParams.variables[paramName][0]
        f_sensParams.close()

    #sensParamValsRow = np.array([[2., 4.]])

    print("\nsensParamValsRow =")
    print(sensParamValsRow)

    # numMetrics x numParams matrix of metric values
    # from sensitivity simulations
    sensMetricValsMatrix = np.zeros((numMetrics, numParams))
    for col in np.arange(numParams):
        f_sens = netcdf.netcdf_file(sensNcFilenames[col], 'r')
        for row in np.arange(numMetrics):
            metricName = metricsNames[row]
            #sensMetricValsMatrix[row,col] = sensMetricValsRefMatrix[row,col]
            sensMetricValsMatrix[row,col] = f_sens.variables[metricName][0]
        f_sens.close()

    #sensMetricValsMatrix = np.array([[1., 2.], [3., 4.]])

    print("\nsensMetricValsMatrix =")
    print(sensMetricValsMatrix)

    f_sensParams.close()

    return(sensParamValsRow, sensMetricValsMatrix)

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
