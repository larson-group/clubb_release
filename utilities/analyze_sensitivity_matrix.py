"""
    Construct and analyze a sensitivity matrix 
    that tells us how user-defined metrics change 
    with changes in parameter values.
"""

def main():

    import numpy as np
    from scipy.io import netcdf
    import pdb
    
    
    paramsNames = np.array(['C5','C8'])
    
    metricsNames = np.array(['SWCF', 'LWCF'])
    
    ncFilenames = np.array(['pathfile1', 'pathfile2'])
    
    if ( len(paramsNames) != len(ncFilenames)   ):
   	print("Number of parameters must equal number of netcdf files.")
   	quit()
    
    # Read netcdf file with metrics and parameters from default simulation
    #f_defaultMetricsParams = netcdf.netcdf_file('canopy/scripts/defaultMetricsParams.nc', 'r')
    # Read netcdf file with changed parameter values from all sensitivity simulations.
    #f_sensParams = netcdf.netcdf_file('canopy/scripts/sensParams.nc', 'r')
    
    # Observed values of our metrics, from, e.g., CERES-EBAF.
    obsMetricValsDict = {'LWCF': 4., 'PRECT': -999., 'SWCF': 6.}
    
    # Values of metrics from the default simulation.
    defaultMetricValsDict = {'SWCF': 5., 'LWCF': 6., 'PRECT': -999.}
    
    # Parameter values from default simulation.
    defaultParamValsDict = {'C15': 42., 'C5': 1., 'C8': 2.}
    
    # Altered parameter values from sensitivity simulations.
    sensParamValsDict = {'C5': 2., 'C15': 42., 'C8': 4.}
    
    # Matrix of metrics from sensitivity simulations.
    # Each column corresponds to one simulation.
    sensMetricValsRefMatrix = np.array([[1., 2.], [3., 4.]])
    
    # Number of tunable parameters
    numParams = len(paramsNames)
    
    print("\nnumParams =")
    print(numParams)
    
    # Number of metrics
    numMetrics = len(metricsNames)
    
    print("\nnumMetrics =")
    print(numMetrics)
    
    # Set up column vector of numMetrics elements containing
    # "true" metric values from observations
    obsMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
   	metricName = metricsNames[idx]
   	obsMetricValsCol[idx] = obsMetricValsDict[metricName]
    
    print("\nobsMetricValsCol =")
    print(obsMetricValsCol)

    defaultMetricValsCol, defaultParamValsRow = \
            setupDefaultVectors(defaultMetricValsDict, metricsNames, 
                                defaultParamValsDict, paramsNames,
                                numParams, numMetrics)
    
    # Create row vector size numParams containing
    # parameter values from sensitivity simulations
    sensParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
   	paramName = paramsNames[idx]
   	sensParamValsRow[0,idx] = sensParamValsDict[paramName]
            # Assume each metric is stored as length-1 array, rather than scalar.  
            #   Hence the "[0]" at the end is needed.
   	#sensParamValsRow[0,idx] = f_sensParams.variables[paramName][0]
    
    #sensParamValsRow = np.array([[2., 4.]])
    
    print("\nsensParamValsRow =")
    print(sensParamValsRow)
    
    # numMetrics x numParams matrix of metric values
    # from sensitivity simulations
    sensMetricValsMatrix = np.zeros((numMetrics, numParams))
    for row in np.arange(numMetrics):
   	metricName = metricsNames[row]
   	for col in np.arange(numParams):
  		ncFilename = ncFilenames[col]
  		sensMetricValsMatrix[row,col] = sensMetricValsRefMatrix[row,col]
    
    #sensMetricValsMatrix = np.array([[1., 2.], [3., 4.]])
    
    print("\nsensMetricValsMatrix =")
    print(sensMetricValsMatrix)
    
    sensMatrix, normlzdSensMatrix = \
        constructSensMatrix(sensParamValsRow, sensMetricValsMatrix, 
                            defaultParamValsRow, defaultMetricValsCol, 
                            obsMetricValsCol, 
                            numParams, numMetrics)
                            
    calcSvd(normlzdSensMatrix)

    #f_defaultMetricsParams.close() 
    #f_sensParams.close()  
        
    print("\nReached the end.")


def constructSensMatrix(sensParamValsRow, sensMetricValsMatrix, 
                        defaultParamValsRow, defaultMetricValsCol, 
                        obsMetricValsCol, 
                        numParams, numMetrics):

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

    import numpy as np                  
                                                      
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
    
    svdInv = np.transpose(vh)*np.diag(sValsTruncInv) *np.transpose(u)
    
    print("\nSVD inverse =")
    print(svdInv)
    
    eigVals, eigVecs = np.linalg.eig(np.transpose(sensMatrix)*sensMatrix)
    
    print("\neigVals =")
    print(eigVals)  
    
    print("\neigVecs = ")
    print(eigVecs)

def setupDefaultVectors(defaultMetricValsDict, metricsNames, 
                        defaultParamValsDict, paramsNames,
                        numParams, numMetrics):    

    import numpy as np            
                                    
    # Set up column vector of numMetrics elements containing
    # metric values from default simulation
    defaultMetricValsCol = np.zeros((numMetrics,1))
    for idx in np.arange(numMetrics):
   	metricName = metricsNames[idx]
   	defaultMetricValsCol[idx] = defaultMetricValsDict[metricName]
            # Assume each metric is stored as length-1 array, rather than scalar.  
            #   Hence the "[0]" at the end is needed.
   	#defaultMetricValsCol[idx] = f_defaultMetricsParams.variables[metricName][0]
    
    print("\ndefaultMetricValsCol =")
    print(defaultMetricValsCol)
    
    # Create row vector size numParams containing
    # parameter values from default simulation
    defaultParamValsRow = np.zeros((1, numParams))
    for idx in np.arange(numParams):
   	paramName = paramsNames[idx]
   	defaultParamValsRow[0,idx] = defaultParamValsDict[paramName]
            # Assume each metric is stored as length-1 array, rather than scalar.  
            #   Hence the "[0]" at the end is needed.
   	#defaultParamValsRow[0,idx] = f_defaultMetricsParams.variables[paramName][0]
    
    #defaultParamValsRow = np.array([[1., 2.]])
    
    print("\ndefaultParamValsRow =")
    print(defaultParamValsRow)
    
    return (defaultMetricValsCol, defaultParamValsRow)

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()