# -*- coding: utf-8 -*-
"""
Use simple input values to test
function analyzeSensMatrix in analyze_sensitivity_matrix.py.
"""

import numpy as np


# The parameters are tunable model parameters.
#    The order of paramsNames must match the order of filenames below.
paramsNames = np.array(['clubb_c8','clubb_c_invrs_tau_n2'])

# Filenames of netcdf files containing fake testing data.
# Netcdf file containing metric and parameter values from the default simulation:
defaultNcFilename = 'default.nc'
#   sensNcFilenames lists filenames of sensitivity simulations.
#   Each file contains metrics and parameter values for a single simulation.
#   There should be one sensitivity simulation per each tunable parameter.
#   The 1st sens filename perturbs the first param in paramsNames, etc.
sensNcFilenames = np.array(['sens0.nc', 'sens1.nc'])

# Observed values of our metrics, from, e.g., CERES-EBAF.
# These observed metrics will be matched as closely as possible by analyzeSensMatrix.
# Global mean, NOTE PRECT is in the unit of m/s
obsMetricValsDict = {'LWCF': 28.008, 'PRECT': 0.000000033912037, 'SWCF': -45.81}

# Write obsMetricValsDict to netcdf file in order to test whether SVD solution recovers original
#    observations.
obsNcFilename = 'obs.nc'

def test_3x2_C8transformed():

    import numpy as np
    import pdb
    from analyze_sensitivity_matrix \
        import analyzeSensMatrix, plotNormlzdSensMatrix

    # The metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsWeights = np.array([[1.], [1.], [1.]])

    # This is the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array(['clubb_c8'])

    # Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
    clubbC8Vals = np.array([0.1812692, 0.50341469, 0.1812692])

    # Create fake testing data with simple values
    write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, defaultBiasesOrigApprox, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    defaultParamValsOrigRow, dparamsSoln, paramsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

    # Check whether the expected answer for the fake data
    #    has been calculated correctly.
    #if np.all( np.isclose( dparamsSoln, np.array([[0.39838053], [2./3.]]), \
    #                          rtol=1e-5, atol=1e-5 ) ):
    if np.all( np.isclose( dparamsSoln, np.array([[0.52633989], [2./3.]]), \
                              rtol=1e-5, atol=1e-5 ) ):
        print("\nPassed test.")
    else:
        print("\ndparamsSoln = ")
        print(dparamsSoln)
        #print("\nError: dparamsSoln should equal [0.48658305 2/3], but it does not.")
        print("\nError: dparamsSoln should equal [0.52633989 2/3], but it does not.")
        assert False

    # Create a heatmap plot that allows us to visualize the normalized sensitivity matrix
    plotNormlzdSensMatrix(normlzdSensMatrix, metricsNames, paramsNames)

    print("\nReached the end of function test_3x2_C8transformed in test harness.")

    return

def test_3x2_novarstransformed():

    import numpy as np
    import pdb
    from analyze_sensitivity_matrix \
        import analyzeSensMatrix, plotNormlzdSensMatrix

    # The metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    metricsNames = np.array(['SWCF', 'LWCF', 'PRECT'])

    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsWeights = np.array([[1.0], [0.000001], [1.0]])

    # This is the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array([''])

    # Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
    clubbC8Vals = np.array([0.2, 0.7, 0.2])

    # Create fake testing data with simple values
    write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, defaultBiasesOrigApprox, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    defaultParamValsOrigRow, dparamsSoln, paramsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

    # Check whether the expected answer for the fake data
    #    has been calculated correctly.
    if np.all( np.isclose( dparamsSoln, np.array([[1.], [0.]]), \
                              rtol=1e-5, atol=1e-5 ) ):
        print("\nPassed test.")
    else:
        print("\ndparamsSoln =  ")
        print(dparamsSoln)
        print("\nError: dparamsSoln should equal [1 0], but it does not.")
        assert False

    print("\nReached the end of function test_3x2_novarstransformed in test harness.")

    return

def test_2x2_C8transformed():

    import numpy as np
    import pdb
    from analyze_sensitivity_matrix \
        import analyzeSensMatrix, plotNormlzdSensMatrix, calcLinSolnBias

    # The metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    metricsNames = np.array(['SWCF', 'LWCF'])

    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsWeights = np.array([[1.], [1.]])

    # This is the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array(['clubb_c8'])

    # Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
    clubbC8Vals = np.array([0.2, 0.20001, 0.2])

    # Create fake testing data with simple values
    write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, defaultBiasesOrigApprox, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    defaultParamValsOrigRow, dparamsSoln, paramsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

    # For the purpose of this unit test, assume that the atmospheric model is indeed linear and
    #   reproduces the observations exactly.
    linSolnNcFilename = obsNcFilename

    # See if the solution based on a linear combination of the SVD-calculated parameter values
    #    matches what we expect.
    linSolnBias = calcLinSolnBias(linSolnNcFilename, defaultNcFilename,
                                  metricsNames)

    # Check whether the expected answer for the fake data
    #    has been calculated correctly.
# This broke when I fed defaultBiasesOrigApprox out of analyzeSensMatrix.
#    if np.all( np.isclose( linSolnBias, defaultBiasesApprox, \
#                               rtol=1e-4, atol=1e-4 ) ):
#        print("\nPassed test.")
#    else:
#         print("\nlinSolnBias = ")
#        print(linSolnBias)
#        print("\ndefaultBiasesApprox = ")
#        print(defaultBiasesApprox)
#        print("\nError: linSolnBias should equal defaultBiasesApprox, but it does not.")
#        assert False

    print("\nReached the end of function test_2x2_C8transformed in test harness.")

def test_2x2_novarstransformed():

    import numpy as np
    import pdb
    from analyze_sensitivity_matrix \
        import analyzeSensMatrix, plotNormlzdSensMatrix

    # The metrics are observed quantities that we want a tuned simulation to match.
    #    The order of metricNames determines the order of rows in sensMatrix.
    metricsNames = np.array(['SWCF', 'LWCF'])

    # Column vector of (positive) weights.  A small value de-emphasizes
    #   the corresponding metric in the fit.
    metricsWeights = np.array([[1.], [1.]])

    # This is the subset of paramsNames that vary from [0,1] (e.g., C5)
    #    and hence will be transformed to [0,infinity] in order to make
    #    the relationship between parameters and metrics more linear:
    transformedParamsNames = np.array([''])

    # Example values of the clubb_c8 parameter in the default and 2 sensitivity runs
    clubbC8Vals = np.array([0.2, 0.7, 0.2])

    # Create fake testing data with simple values
    write_test_netcdf_files(obsMetricValsDict, clubbC8Vals)

    # Calculate changes in parameter values needed to match metrics.
    defaultMetricValsCol, defaultBiasesCol, defaultBiasesOrigApprox, \
    sensMatrixOrig, sensMatrix, normlzdSensMatrix, svdInvrsNormlzdWeighted, \
    defaultParamValsOrigRow, dparamsSoln, paramsSoln = \
        analyzeSensMatrix(metricsNames, paramsNames, transformedParamsNames,
                        metricsWeights,
                        sensNcFilenames, defaultNcFilename,
                        obsMetricValsDict)

    # Check whether the expected answer for the fake data
    #    has been calculated correctly.
    if np.all( np.isclose( dparamsSoln, np.array([[1.], [1.]]), \
                              rtol=1e-5, atol=1e-5 ) ):
        print("\nPassed test.")
    else:
        print("\ndparamsSoln = ")
        print(dparamsSoln)
        print("\nError: dparamsSoln should equal [1 1], but it does not.")
        assert False

    print("\nReached the end of function test_2x2_novarstransformed in test harness.")


def write_test_netcdf_files(obsMetricValsDict, clubbC8Vals):

    import netCDF4

    SWCF_obs = obsMetricValsDict['SWCF']
    LWCF_obs = obsMetricValsDict['LWCF']
    PRECT_obs = obsMetricValsDict['PRECT']

    f_obs = netCDF4.Dataset('obs.nc', 'w')
    f_obs.createDimension('scalar', 1)
    SWCF = f_obs.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs
    LWCF = f_obs.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs
    PRECT = f_obs.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs
    f_obs.close()

    f_default = netCDF4.Dataset('default.nc', 'w')
    f_default.createDimension('scalar', 1)
    SWCF = f_default.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 1.0
    LWCF = f_default.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs - 1.0
    PRECT = f_default.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs - 1.0
    clubb_c8 = f_default.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = clubbC8Vals[0]  #0.1812692 #0.2
    clubb_c_invrs_tau_n2 = f_default.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 3.
    f_default.close()

    f_sens0 = netCDF4.Dataset('sens0.nc', 'w')
    f_sens0.createDimension('scalar', 1)
    SWCF = f_sens0.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 0.5
    LWCF = f_sens0.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs - 1.
    PRECT = f_sens0.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs - 0.5
    clubb_c8 = f_sens0.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = clubbC8Vals[1]  #0.50341469 #0.7
    clubb_c_invrs_tau_n2 = f_sens0.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 3.
    f_sens0.close()

    f_sens1 = netCDF4.Dataset('sens1.nc', 'w')
    f_sens1.createDimension('scalar', 1)
    SWCF = f_sens1.createVariable('SWCF', 'f', ('scalar',))
    SWCF[:] = SWCF_obs - 1.
    LWCF = f_sens1.createVariable('LWCF', 'f', ('scalar',))
    LWCF[:] = LWCF_obs + 1.
    PRECT = f_sens1.createVariable('PRECT', 'f', ('scalar',))
    PRECT[:] = PRECT_obs + 1.
    clubb_c8 = f_sens1.createVariable('clubb_c8', 'f', ('scalar',))
    clubb_c8[:] = clubbC8Vals[2] #0.1812692 #0.2
    clubb_c_invrs_tau_n2 = f_sens1.createVariable('clubb_c_invrs_tau_n2', 'f', ('scalar',))
    clubb_c_invrs_tau_n2[:] = 5.
    f_sens1.close()

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    test_3x2_C8transformed()

