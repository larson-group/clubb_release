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
        'sens0.nc', \
        'sens1.nc'
              ])
    sensNcFilenamesExt = \
    np.array([ \
        'sens0.nc', \
        'sens1.nc'
              ])

    # Observed values of our metrics, from, e.g., CERES-EBAF.
    # These observed metrics will be matched as closely as possible by analyzeSensMatrix.
    # Global mean, NOTE PRECT is in the unit of m/s
    obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, \
    defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC, \
    normlzdWeightedDefaultBiasesApprox, normlzdWeightedDefaultBiasesApproxPC, \
    defaultBiasesOrigApprox, defaultBiasesOrigApproxPC, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, \
    normlzdWeightedSensMatrix, biasNormlzdSensMatrix, svdInvrsNormlzdWeighted, \
    vhNormlzd, uNormlzd, sNormlzd, \
    magParamValsRow, \
    defaultParamValsOrigRow, dparamsSoln, dnormlzdParamsSoln, \
    dparamsSolnPC, dnormlzdParamsSolnPC, \
    paramsSoln, paramsLowVals, paramsHiVals, \
    paramsSolnPC, paramsLowValsPC, paramsHiValsPC = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                      metricsWeights,
                      sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
                       obsMetricValsDict)

    # See if new global simulation output based on a linear combination
    #    of the SVD-calculated parameter values matches what we expect.
    linSolnDiff = calcLinSolnDiff(linSolnNcFilename, defaultNcFilename,
                                  metricsNames)

    # Create a heatmap plot that allows us to visualize the normalized sensitivity matrix
    plotNormlzdSensMatrix(normlzdSensMatrix, metricsNames, paramsNames)

    print("\nReached the end of main.")

    return

def analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                      metricsWeights,
                      sensNcFilenames, sensNcFilenamesExt, defaultNcFilename,
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
    from sklearn.preprocessing import normalize

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
            setupObsCol(obsMetricValsDict, metricsNames)

    print("\nobsMetricValsCol =")
    print(obsMetricValsCol)

    # Set up a column vector of metric values from the default simulation
    defaultMetricValsCol = \
        setupDefaultMetricValsCol(metricsNames, defaultNcFilename)

    print("\ndefaultMetricValsCol =")
    print(defaultMetricValsCol)

#    sys.exit("Prints")

    # Based on the default simulation,
    #    set up a column vector of metrics and a row vector of parameter values.
    defaultParamValsRow, defaultParamValsOrigRow = \
            setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
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
                print("\nSensitivity filename =", sensNcFilenames[sensFileIdx])
                print("sensFileIdx =", sensFileIdx)
                print("Parameter name = ", paramName)
                print("paramIdx =", paramIdx, ", defaultParamVal =", defaultParamVal, ", sensParamVal =", sensParamVal)
                sys.exit("Error: a sensitivity simulation has left a designated parameter at its default value.")
            if sensFileIdx!=paramIdx and not sensDefaultAreClose:
                print("\nSensitivity filename =", sensNcFilenames[sensFileIdx])
                print("sensFileIdx =", sensFileIdx)
                print("Parameter name = ", paramName)
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
    sensMetricValsMatrixExt, sensParamValsRowExt, sensParamValsOrigRowExt = \
            setupSensArrays(metricsNames, paramsNames, transformedParamsNames,
                    numMetrics, numParams,
                    sensNcFilenamesExt)

    #pdb.set_trace()

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
    
    # Matrix of metric values from default simulation
    # Each column in the matrix is repeated numParams times, for later multiplication
    defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    # Use untransformed (original) parameter values.
    defaultBiasesColOrig, sensMatrixOrig, normlzdSensMatrixOrig, biasNormlzdSensMatrixOrig = \
         constructSensMatrix(sensMetricValsMatrix, sensParamValsOrigRow,
                            sensMetricValsMatrixExt, sensParamValsOrigRowExt, # deriv = sensExt - sens
                            #defaultMetricValsMatrix, defaultParamValsOrigRow, # deriv = sens - default
                            defaultMetricValsCol,
                            np.full_like(magParamValsRow,1.0),
                            obsMetricValsCol,
                            numMetrics, numParams,
                            beVerbose=False)

    print("\nsensMatrixOrig =")
    print(sensMatrixOrig)

    # Calculate the sensitivity matrix and the sensitivity matrix
    # normalized by the discrepancies from observations in default simulation.
    # Use transformed parameter values.
    defaultBiasesCol, sensMatrix, normlzdSensMatrix, biasNormlzdSensMatrix = \
         constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                            sensMetricValsMatrixExt, sensParamValsRowExt, # deriv = sensExt - sens
                            #defaultMetricValsMatrix, defaultParamValsRow, # deriv = sens - default
                            defaultMetricValsCol,
                            magParamValsRow,
                            obsMetricValsCol,
                            numMetrics, numParams,
                            beVerbose=True)

    uNormlzd, sNormlzd, vhNormlzd = np.linalg.svd(normlzdSensMatrix, full_matrices=False)

    # In order to weight certain metrics, multiply each row of normlzdSensMatrix
    # by metricsWeights
    normlzdWeightedSensMatrix = np.diag(np.transpose(metricsWeights)[0]) @ normlzdSensMatrix

    # Compute dot products of  bias vector and matrix columns
    normlzdDefaultBiasesCol = ( metricsWeights * (-defaultBiasesCol) /
                                np.abs(obsMetricValsCol) )
    augSensMatrix = np.hstack( (normlzdDefaultBiasesCol, normlzdWeightedSensMatrix) )
    augSensMatrixColNorm = normalize(augSensMatrix, axis=0, norm='l2')
    # Normalized this way, dotMatrix = cos(angle between columns).
    # The first column is the normalized bias; the other columns are the sensitivity matrix.
    dotMatrix = augSensMatrixColNorm.T @ augSensMatrixColNorm

    print("\ndotMatrix =")
    print(dotMatrix)

    # Calculate lower bound on normalized parameter perturbations
    sensMatrixRowMag = np.linalg.norm(normlzdWeightedSensMatrix, axis=1)
    dpMin = np.abs(normlzdDefaultBiasesCol) / np.atleast_2d(sensMatrixRowMag).T

    print("\nsensMatrixRowMag =")
    print(np.atleast_2d(sensMatrixRowMag).T)

    print("\nnormlzdDefaultBiasesCol =")
    print(normlzdDefaultBiasesCol)

    print("\n U*b =")
    print(np.atleast_2d(sensMatrixRowMag).T * normlzdDefaultBiasesCol)

    print("\ndpMin =")
    print(dpMin)
    #pdb.set_trace()

    # sValsRatio = a threshold ratio of largest singular value
    #              to the smallest retained singular value.
    # If sValsRatio is large enough, then all singular vectors will be kept.
    # If sValsRatio is 1, then only the first singular vector will be kept.
    sValsRatio = 30.

    # Calculate inverse of the singular value decomposition.
    # This gives the recommended changes to parameter values.
    svdInvrsNormlzdWeighted, svdInvrsNormlzdWeightedPC, \
    sValsTruncInvNormlzdWeighted, sValsTruncInvNormlzdWeightedPC, \
    vhNormlzdWeighted, uNormlzdWeighted, sNormlzdWeighted = \
    calcSvdInvrs(normlzdWeightedSensMatrix, sValsRatio)

    #print("\nNormalized, weighted SVD inverse =")
    #print(svdInvrsNormlzdWeighted)

    paramsSoln, paramsLowVals, paramsHiVals, dparamsSoln, dnormlzdParamsSoln, \
    defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals = \
    calcParamsSoln(svdInvrsNormlzdWeighted, metricsWeights, magParamValsRow, \
                   sensMatrix, normlzdWeightedSensMatrix, \
                   obsMetricValsCol, defaultBiasesCol, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeighted, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames )

    paramsSolnPC, paramsLowValsPC, paramsHiValsPC, dparamsSolnPC, dnormlzdParamsSolnPC, \
    defaultBiasesApproxPC, defaultBiasesApproxLowValsPC, defaultBiasesApproxHiValsPC = \
    calcParamsSoln(svdInvrsNormlzdWeightedPC, metricsWeights, magParamValsRow, \
                   sensMatrix, normlzdWeightedSensMatrix, \
                   obsMetricValsCol, defaultBiasesCol, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeightedPC, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames )

    print("\ndefaultBiasesApprox =")
    print(defaultBiasesApprox)

    print("\nmag(dnormlzdParamsSoln) =")
    print(np.linalg.norm(dnormlzdParamsSoln))

    print("\nmag(dnormlzdParamsSolnPC) =")
    print(np.linalg.norm(dnormlzdParamsSolnPC))

    print("\ndparamsSoln =")
    print(dparamsSoln)

    print("\nparamsSoln =")
    print(paramsSoln)

    # This check is currently broken if any params span [0,1], e.g., C5
    #if (transformedParamsNames == np.array([''])).all():
    print("\nCheck: Does defaultBiasesOrigApprox approximate -defaultBiasesCol above?")
    defaultBiasesOrigApprox = sensMatrixOrig @ dparamsSoln
    print("defaultBiasesOrigApprox =")
    print(defaultBiasesOrigApprox)

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    normlzdWeightedDefaultBiasesApprox = \
            normlzdWeightedSensMatrix @ ( dparamsSoln / np.transpose(magParamValsRow) )

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    normlzdWeightedDefaultBiasesApproxPC = \
            normlzdWeightedSensMatrix @ ( dparamsSolnPC / np.transpose(magParamValsRow) )

    defaultBiasesOrigApproxPC = sensMatrixOrig @ dparamsSolnPC

    #pdb.set_trace()

    return (defaultMetricValsCol, defaultBiasesCol, \
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
            paramsSolnPC, paramsLowValsPC, paramsHiValsPC)

def constructSensMatrix(sensMetricValsMatrix, sensParamValsRow,
                        defaultMetricValsMatrix, defaultParamValsRow,
                        defaultMetricValsCol,
                        magParamValsRow,
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
    #defaultMetricValsMatrix = defaultMetricValsCol @ np.ones((1,numParams))

    #print("\ndefaultMetricValsMatrix =")
    #print(defaultMetricValsMatrix)

    # Sensitivity simulation metrics minus default simulation metrics
    #pdb.set_trace()
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

    # Sensitivity matrix of derivatives, dmetrics/dparams.
    # sensMatrix is not normalized or weighted by any factor.
    sensMatrix = dmetricsMatrix * invrsDparamsMatrix

    if beVerbose:
        print("\nsensMatrix =")
        print(sensMatrix)

    # Store biases in default simulation
    # defaultBiasesCol = + delta_b
    defaultBiasesCol = np.subtract(defaultMetricValsCol, obsMetricValsCol)


    if beVerbose:
        print("\ndefaultBiasesCol =")
        print(defaultBiasesCol)

    # Matrix of inverse biases.
    # Used for forming normalized sensitivity derivatives.
    #pdb.set_trace()
    invrsObsMatrix = np.reciprocal(np.abs(obsMetricValsCol)) @ np.ones((1,numParams))

#    if beVerbose:
#        print("\ninvrsBiasesMatrix =")
#        print(invrsObsMatrix)

    # Form matrix of default parameter values, for later normalization of the sensitivity matrix
    magParamValsMatrix = np.ones((numMetrics,1)) @ magParamValsRow

#    if beVerbose:
#        print("\nmagParamValsMatrix =")
#        print(magParamValsMatrix)

    # Sensitivity matrix, normalized by obs and parameter values, but not weighted
    normlzdSensMatrix = sensMatrix * invrsObsMatrix * magParamValsMatrix

    # Matrix of inverse biases.
    # Used for forming sensitivity matrix normalized by biases, not observed values.
    invrsBiasesMatrix = np.reciprocal(np.abs(defaultBiasesCol)) @ np.ones((1,numParams))

    # Sensitivity matrix, normalized by obs and parameter values, but not weighted
    biasNormlzdSensMatrix = sensMatrix * invrsBiasesMatrix * magParamValsMatrix

    if beVerbose:
        print("\nnormlzdSensMatrix =")
        print(normlzdSensMatrix)

    return  (defaultBiasesCol, sensMatrix, normlzdSensMatrix, biasNormlzdSensMatrix)


def calcSvdInvrs(normlzdWeightedSensMatrix, sValsRatio):
    """
    Input: sensitivity matrix
    Output: singular value decomposition of sensitivity matrix
        and related quantities, like the inverse of the SVD.
    """


    import numpy as np
    import sys
    import pdb

    # vh = V^T = transpose of right-singular vector matrix, V.
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

    # Assertion check: is svdInvrs truly the inverse of normlzdWeightedSensMatrix?
    numParams = normlzdWeightedSensMatrix.shape[1] # = number of columns
    if not np.all( np.isclose(np.identity(numParams), svdInvrs @ normlzdWeightedSensMatrix, \
                      rtol=2e-6 , atol=2e-6 ) ):
        print("\nsvdInvrs @ normlzdWeightedSensMatrix =")
        print(svdInvrs @ normlzdWeightedSensMatrix)
        sys.exit("\nError: svdInvrs is not the inverse of normlzdWeightedSensMatrix")

    #print("\nSVD inverse =")
    #print(svdInvrs)


    # Delete the small singular values in order to show just the most important patterns.
    # After this deletion, store inverse singular values in sValsInvPC.
    # Syntax: (a_truncated = a[0:a.size-2] lops off the last two elements of `a`.)
    # sValsInvPC[0:sValsTruncInv.size-2] = sValsTruncInv[0:sValsTruncInv.size-2]
    sValsInvPC = np.zeros_like(sValsTruncInv)
    for idx, val in np.ndenumerate(sValsTruncInv):
        # If a singular value is much smaller than largest singular value, then zero it out.
        if val/sValsTruncInv[0] > sValsRatio:
            sValsInvPC[idx] = 0.
        else:
            sValsInvPC[idx] = sValsTruncInv[idx]

    print("\nsValsInvPC =")
    print(sValsInvPC)

    #pdb.set_trace()

    svdInvrsPC = np.transpose(vh) @ np.diag(sValsInvPC) @ np.transpose(u)

    # eigVals, eigVecs = np.linalg.eig(np.transpose(normlzdWeightedSensMatrix) @ normlzdWeightedSensMatrix)

    #print("\neigVals =")
    #print(eigVals)

    #print("\neigVecs = ")
    #print(eigVecs)

    return ( svdInvrs, svdInvrsPC, sValsTruncInv, sValsInvPC, vh, u, s )


def calcParamsSoln(svdInvrsNormlzdWeighted, metricsWeights, magParamValsRow, \
                   sensMatrix, normlzdWeightedSensMatrix, \
                   obsMetricValsCol, defaultBiasesCol, defaultParamValsOrigRow, \
                   sValsTruncInvNormlzdWeighted, vhNormlzdWeighted, \
                   numParams, paramsNames, transformedParamsNames ):

    import numpy as np
    import pdb
    import sys

    # Calculate solution using transformed SVD matrix.  However, dparamsSoln is not normalized.
    # dnormlzdParamsSoln = dparamsSoln / ( normalizing param value from default simulation
    #                                       (or sens sim if default param value is zero.) )
    dnormlzdParamsSoln = svdInvrsNormlzdWeighted @ \
                      ( metricsWeights * (-defaultBiasesCol) / np.abs(obsMetricValsCol) )
    # Unnormalized parameter value perturbations from the default values, delta_p
    dparamsSoln = dnormlzdParamsSoln * np.transpose(magParamValsRow)
    # defaultBiasesApprox = (forward model soln - default soln)
    defaultBiasesApprox = ( normlzdWeightedSensMatrix @ dnormlzdParamsSoln ) \
                          * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol)

    # As a check, do an alternative calculation of defaultBiasesApprox based on
    #    the unnormalized, unweighted sensitivity matrix:
    defaultBiasesApproxAlt = sensMatrix @ dparamsSoln
    if ( not np.all( np.isclose(defaultBiasesApproxAlt, defaultBiasesApprox) ) ):
        print("\ndefaultBiasesApprox as computed from sensMatrix =")
        print(defaultBiasesApproxAlt)
        print("\ndefaultBiasesApprox as computed from normlzdWeightedSensMatrix =")
        print(defaultBiasesApprox)
        sys.exit("\nError: Two calculations of defaultBiasesApprox do not agree.")

    #pdb.set_trace()
    paramsSoln = np.transpose(defaultParamValsOrigRow) + dparamsSoln
    # Create matrix whose columns are the columns of v divided by singular values
    #    (see 15.4.18 of Numerical Recipes)
    dparamsErrMatrix = np.transpose(vhNormlzdWeighted) * sValsTruncInvNormlzdWeighted
    # Create a column vector which contains the standard deviation of error for each parameter
    #    difference
    dparamsErrStdNormlzd = np.sqrt( np.reshape( np.sum( dparamsErrMatrix * dparamsErrMatrix, axis=1 ), (-1, 1)) )
    # Rescale errors by parameter magnitude scales
    dparamsErrStd = dparamsErrStdNormlzd * np.transpose(magParamValsRow)
    # Do an alternative calculation for the error bars based on normal equations.
    # Follow Section 15.4.1 of Numerical Recipes
    covar_matrix = np.linalg.inv( normlzdWeightedSensMatrix.T @ normlzdWeightedSensMatrix )
    dparamsErrStdNormlzdCov = np.sqrt( covar_matrix.diagonal() )
    dparamsErrStdCov = dparamsErrStdNormlzdCov * np.transpose(magParamValsRow)

    dparamsLowVals = dparamsSoln - dparamsErrStd
    paramsLowVals = np.transpose(defaultParamValsOrigRow) + dparamsLowVals
    dparamsHiVals = dparamsSoln + dparamsErrStd
    paramsHiVals = np.transpose(defaultParamValsOrigRow) + dparamsHiVals
    defaultBiasesApproxLowVals = sensMatrix @ dparamsLowVals
    defaultBiasesApproxHiVals = sensMatrix @ dparamsHiVals

    #pdb.set_trace()
    # Transform some variables from [-inf,inf] back to [0,inf] range
    for idx in np.arange(numParams):
        paramName = paramsNames[idx]
        if paramName in transformedParamsNames:
            #paramsSoln[idx,0] = 1.0-np.exp(-paramsSoln[idx,0])
            # The equation below means:
            # param_transformed_updated = delta_param_transformed + param_transformed_orig
            #    where param_transformed_orig = log( param_untransformed_orig )
            paramsSoln[idx,0] = np.exp(dparamsSoln[idx,0]) * defaultParamValsOrigRow[0,idx]
            dparamsSoln[idx,0] = paramsSoln[idx,0] - defaultParamValsOrigRow[0,idx]
            paramsLowVals[idx,0] = np.exp(dparamsLowVals[idx,0]) * defaultParamValsOrigRow[0,idx]
            paramsHiVals[idx,0] = np.exp(dparamsHiVals[idx,0]) * defaultParamValsOrigRow[0,idx]

    return (paramsSoln, paramsLowVals, paramsHiVals, dparamsSoln, dnormlzdParamsSoln, \
            defaultBiasesApprox, defaultBiasesApproxLowVals, defaultBiasesApproxHiVals)

def setupObsCol(obsMetricValsDict, metricsNames):
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

    print("\nobsMetricValsCol =")
    print(obsMetricValsCol)

    return obsMetricValsCol

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

    print("\ndefaultMetricValsCol =")
    print(defaultMetricValsCol)

    f_defaultMetricsParams.close()

    return defaultMetricValsCol

def setupDefaultParamVectors(metricsNames, paramsNames, transformedParamsNames,
                        numMetrics, numParams,
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

    print("\ndefaultParamValsOrigRow =")
    print(defaultParamValsOrigRow)

    print("\ndefaultParamValsRow =")
    print(defaultParamValsRow)

    f_defaultMetricsParams.close()

    return (defaultParamValsRow, defaultParamValsOrigRow)

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


def calcLinSolnDiff(linSolnNcFilename, defaultNcFilename,
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

    linSolnDiff = linSolnMetricValsCol - defaultMetricValsCol

    print("\nlinSolnDiff")
    print(linSolnDiff)

    return linSolnDiff

def findOutliers(normlzdSensMatrix, normlzdWeightedSensMatrix, \
                 defaultBiasesCol, obsMetricValsCol, magParamValsRow, defaultParamValsOrigRow):

    import numpy as np
    from sklearn import linear_model
    import pdb

#    ransac = linear_model.RANSACRegressor(max_trials=1000,random_state=0,
#                                          base_estimator=linear_model.LinearRegression(), residual_threshold=None)
#    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) )
#    defaultBiasesApproxRansac = ransac.predict(normlzdSensMatrix) * np.abs(obsMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( ransac.estimator_.coef_ )
#    ransac = linear_model.RANSACRegressor(max_trials=10000,random_state=0,
#             base_estimator=linear_model.ElasticNet(fit_intercept=False, random_state=0, tol=1e-3,
#                                                    l1_ratio=0.0, alpha=5),
#                                          residual_threshold=None)
#    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) )
#    defaultBiasesApproxRansac = np.transpose( np.atleast_2d( ransac.predict(normlzdSensMatrix) ) ) \
#                                * np.abs(obsMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( np.atleast_2d( ransac.estimator_.coef_ ) )
#    inlier_mask = ransac.inlier_mask_
#    outlier_mask = np.logical_not(inlier_mask)


    ransac = linear_model.HuberRegressor(fit_intercept=False)
    ransac.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) )
    defaultBiasesApproxRansac = np.transpose( np.atleast_2d( ransac.predict(normlzdSensMatrix) ) ) \
                                * np.abs(obsMetricValsCol)
    dnormlzdParamsSolnRansac = np.transpose( np.atleast_2d( ransac.coef_ ) )
#    defaultBiasesApproxRansac = ransac.predict(normlzdSensMatrix) * np.abs(obsMetricValsCol)
#    dnormlzdParamsSolnRansac = np.transpose( ransac.coef_ )

    dparamsSolnRansac = dnormlzdParamsSolnRansac * np.transpose(magParamValsRow)
    paramsSolnRansac = np.transpose(defaultParamValsOrigRow) + dparamsSolnRansac

    outlier_mask = ransac.outliers_
    inlier_mask = np.logical_not(outlier_mask)

    print( "paramsSolnRansac = ", paramsSolnRansac )
    print( "dparamsSolnRansac = ", dparamsSolnRansac )

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    normlzdWeightedDefaultBiasesApproxRansac = \
            normlzdWeightedSensMatrix @ dnormlzdParamsSolnRansac

    #pdb.set_trace()

    return (outlier_mask, defaultBiasesApproxRansac, normlzdWeightedDefaultBiasesApproxRansac, \
            dnormlzdParamsSolnRansac, paramsSolnRansac)

def findParamsUsingElastic(normlzdSensMatrix, normlzdWeightedSensMatrix, \
                 defaultBiasesCol, obsMetricValsCol, metricsWeights, \
                 magParamValsRow, defaultParamValsOrigRow, \
                 normlzdCurvMatrix):

    import numpy as np
    from sklearn.linear_model import ElasticNet, ElasticNetCV
    from sklearn.linear_model import Lasso, LassoCV
    import pdb

    #regr = ElasticNet(fit_intercept=True, random_state=0, tol=1e-10, l1_ratio=0.5, alpha=0.01)
    #regr =Lasso(fit_intercept=True, random_state=0, tol=1e-10, alpha=0.01) # don't fit intercept!;use line below
    regr = Lasso(fit_intercept=False, random_state=0, tol=1e-10, alpha=0.01)
    #regr = LassoCV(fit_intercept=True, random_state=0, eps=1e-5, tol=1e-10, cv=metricsWeights.size)
    #print( "alpha_ = ", regr.alpha_ )
    regr.fit(normlzdWeightedSensMatrix, -metricsWeights * defaultBiasesCol / np.abs(obsMetricValsCol) )
    #regr.fit(normlzdSensMatrix, -defaultBiasesCol / np.abs(obsMetricValsCol) )

    defaultBiasesApproxElastic = np.transpose( np.atleast_2d(regr.predict(normlzdWeightedSensMatrix)) ) \
                          * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol)
    dnormlzdParamsSolnElastic = np.transpose( np.atleast_2d(regr.coef_) )
    dparamsSolnElastic = dnormlzdParamsSolnElastic * np.transpose(magParamValsRow)
    paramsSolnElastic = np.transpose(defaultParamValsOrigRow) + dparamsSolnElastic

    print( "paramsSolnElastic = ", paramsSolnElastic )
    print( "dparamsSolnElastic = ", dparamsSolnElastic )

    #pdb.set_trace()

    # If the solution were perfect, this variable would equal
    #     the normalized, weighted right-hand side.
    defaultBiasesApproxElasticNonlin = \
            normlzdWeightedSensMatrix @ dnormlzdParamsSolnElastic \
                        * np.reciprocal(metricsWeights) * np.abs(obsMetricValsCol) \
            + 0.5 * normlzdCurvMatrix @ (dnormlzdParamsSolnElastic**2) * np.abs(obsMetricValsCol)

    #pdb.set_trace()

    return (defaultBiasesApproxElastic, defaultBiasesApproxElasticNonlin, \
            dnormlzdParamsSolnElastic, paramsSolnElastic)

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
